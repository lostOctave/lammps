// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing authors: Eugen Rozic (University College London)
------------------------------------------------------------------------- */

#include "pair_wca_att.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include <iostream>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;
/* ---------------------------------------------------------------------- */

PairWCAatt::PairWCAatt(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairWCAatt::~PairWCAatt()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon);
    memory->destroy(eps_wca);
    memory->destroy(sigma);
    memory->destroy(w);
    memory->destroy(cut);
    memory->destroy(wcaflag);

    memory->destroy(lj12_e);
    memory->destroy(lj6_e);
    memory->destroy(lj12_f);
    memory->destroy(lj6_f);
    memory->destroy(local_force_x);
    memory->destroy(local_force_y);
    memory->destroy(force_x);
    memory->destroy(force_y);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairWCAatt::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  int n_atoms = atom->natoms;
  memory->create(setflag, n+1, n+1, "pair:setflag");
  memory->create(cutsq, n+1, n+1, "pair:cutsq");

  memory->create(cut, n+1, n+1, "pair:cut");
  memory->create(epsilon, n+1, n+1, "pair:epsilon");
  memory->create(eps_wca, n+1, n+1, "pair:epsilon");
  memory->create(sigma, n+1, n+1, "pair:sigma");
  memory->create(w, n+1, n+1, "pair:w");
  memory->create(wcaflag, n+1, n+1, "pair:wcaflag");

  memory->create(lj12_e, n+1, n+1, "pair:lj12_e");
  memory->create(lj6_e, n+1, n+1, "pair:lj6_e");
  memory->create(lj12_f, n+1, n+1, "pair:lj12_f");
  memory->create(lj6_f, n+1, n+1, "pair:lj6_f");
  memory->create(local_force_x, n_atoms+1,"pair:local_force_x");
  memory->create(local_force_y, n_atoms+1,"pair:local_force_y");
  memory->create(force_x, n_atoms+1,"pair:force_x");
  memory->create(force_y, n_atoms+1,"pair:force_y");

  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
      wcaflag[i][j] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairWCAatt::settings(int narg, char **arg)
{
  if (narg != 1) {
    error->all(FLERR, "Illegal pair_style command (wrong number of params)");
  }

  cut_global = utils::numeric(FLERR, arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut[i][j] = cut_global;
  }
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairWCAatt::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 6)
    error->all(FLERR, "Incorrect args for pair coefficients (too few or too many)");

  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR, arg[3],false,lmp);
  double cut_one = utils::numeric(FLERR, arg[4],false,lmp);
  double epsilon_wca = utils::numeric(FLERR, arg[5],false,lmp);

  double wca_one = 0;
  //int flag;
  //char cheps_atom[] = "eps_atom";
  //int index_eps_atom = atom->find_custom(cheps_atom,flag,1);
  //double *eps_atom = atom->dvector[index_eps_atom];

  if (cut_one < sigma_one) {
    error->all(FLERR, "Incorrect args for pair coefficients (cutoff < sigma)");
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      eps_wca[i][j] = epsilon_wca;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      wcaflag[i][j] = wca_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR, "Incorrect args for pair coefficients (none set)");
}

/* ----------------------------------------------------------------------
   init specific to this pair style (unnecessary)
------------------------------------------------------------------------- */

/*
void PairWCAatt::init_style()
{
  neighbor->request(this,instance_me);
}
*/

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairWCAatt::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, "Mixing not supported in pair_style cosine/squared");

  epsilon[j][i] = epsilon[i][j];
  eps_wca[j][i] = eps_wca[i][j];
  sigma[j][i] = sigma[i][j];
  cut[j][i] = cut[i][j];

  w[j][i] = w[i][j] = cut[i][j] - sigma[i][j];

  lj12_e[j][i] = lj12_e[i][j] = eps_wca[i][j] * pow(sigma[i][j], 12.0);
  lj6_e[j][i] = lj6_e[i][j] = 2.0 * eps_wca[i][j] * pow(sigma[i][j], 6.0);
  lj12_f[j][i] = lj12_f[i][j] = 12.0 * eps_wca[i][j] * pow(sigma[i][j], 12.0);
  lj6_f[j][i] = lj6_f[i][j] = 12.0 * eps_wca[i][j] * pow(sigma[i][j], 6.0);

  // Note: cutsq is set in pair.cpp

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   this is here to throw errors & warnings for given options
------------------------------------------------------------------------- */

void PairWCAatt::modify_params(int narg, char **arg)
{
  Pair::modify_params(narg, arg);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mix") == 0) {
      error->all(FLERR, "pair_modify mix not supported for pair_style cosine/squared");
    } else if (strcmp(arg[iarg], "shift") == 0) {
      error->warning(FLERR, "pair_modify shift has no effect on pair_style cosine/squared");
      offset_flag = 0;
    } else if (strcmp(arg[iarg], "tail") == 0) {
      error->warning(FLERR, "pair_modify tail has no effect on pair_style cosine/squared");
      tail_flag = 0;
    }
    iarg++;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWCAatt::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&eps_wca[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWCAatt::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0)
        utils::sfread(FLERR,&setflag[i][j], sizeof(int), 1, fp,nullptr,error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon[i][j], sizeof(double), 1, fp,nullptr,error);
          utils::sfread(FLERR,&eps_wca[i][j], sizeof(double), 1, fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j], sizeof(double), 1, fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j], sizeof(double), 1, fp,nullptr,error);
        }
	MPI_Bcast(&eps_wca[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWCAatt::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWCAatt::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cut_global, sizeof(double), 1, fp,nullptr,error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairWCAatt::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %g\n", i, epsilon[i][i], sigma[i][i], cut[i][i], eps_wca[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairWCAatt::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g\n", i, j, epsilon[i][j], sigma[i][j],
              cut[i][j], eps_wca[i][j] );
}

/* ---------------------------------------------------------------------- */

void PairWCAatt::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double r, rsq, r2inv, r6inv;
  double factor_lj, force_lj, force_cos, cosone;

  ev_init(eflag, vflag);
  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int n_atoms = atom->natoms;
  int *atom_id = atom->tag;
  //cout<<n_atoms<<endl;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  int flag;
  int cols;

  char cheps_atom[] = "eps_atom";
  int index_eps_atom = atom->find_custom(cheps_atom,flag,cols);
  double *eps_atom = atom->dvector[index_eps_atom];
  // here we store pair wise interaction forces on atoms
  char chforce_atom_x[] = "force_pair_x";
  int index_force_atom_x = atom->find_custom(chforce_atom_x,flag,cols);
  double *force_atom_x = atom->dvector[index_force_atom_x];

  char chforce_atom_y[] = "force_pair_y";
  int index_force_atom_y = atom->find_custom(chforce_atom_y,flag,cols);
  double *force_atom_y = atom->dvector[index_force_atom_y];

  double epsilon_att_pair;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  for (int i_all =0; i_all<n_atoms;i_all++) {
    local_force_x[i_all] = 0;
    local_force_y[i_all] = 0;
    force_x[i_all] = 0;
    force_y[i_all] = 0;
  }
  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    //cout<<i<<" "<<atom_id[i]<<endl;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      //epsilon_att_pair = sqrt(eps_atom[i]*eps_atom[j]);
      if (rsq < cutsq[itype][jtype]) {

        /*
          This is exactly what the "single" method does, in fact it could be called
          here instead of repeating the code but here energy calculation is optional
          so a little bit of calculation is possibly saved
        */

        r = sqrt(rsq);
	//cout<<i<<"	"<<j<<"	"<<r<<"	"<<epsilon[itype][jtype]<<endl;

        if (r <= sigma[itype][jtype]) {
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          force_lj = r6inv*(lj12_f[itype][jtype]*r6inv - lj6_f[itype][jtype]);
          fpair = factor_lj*force_lj*r2inv;
          evdwl = factor_lj*r6inv*(lj12_e[itype][jtype]*r6inv - lj6_e[itype][jtype]);
          //evdwl += factor_lj*(eps_wca[itype][jtype]-epsilon_att_pair); // energy is zero at cutoff
          evdwl += factor_lj*(eps_wca[itype][jtype]-epsilon[itype][jtype]); // energy is zero at cutoff
        }else{
          force_cos = -(MY_PI*epsilon[itype][jtype] / (2.0*w[itype][jtype])) *
                      sin(MY_PI*(r-sigma[itype][jtype]) / w[itype][jtype]);
          fpair = factor_lj*force_cos / r;
          cosone = cos(MY_PI*(r-sigma[itype][jtype]) / (2.0*w[itype][jtype]));
          evdwl = -factor_lj*epsilon[itype][jtype]*cosone*cosone;
	      /*force_cos = -(MY_PI*epsilon_att_pair / (2.0*w[itype][jtype])) *
                      sin(MY_PI*(r-sigma[itype][jtype]) / w[itype][jtype]);
          fpair = factor_lj*force_cos / r;
          cosone = cos(MY_PI*(r-sigma[itype][jtype]) / (2.0*w[itype][jtype]));
          evdwl = -factor_lj*epsilon_att_pair*cosone*cosone;*/
	}
	//cout << i <<" "<<j<<" "<<endl;
	local_force_x[atom_id[i]-1] += delx*fpair;
	local_force_y[atom_id[i]-1] += dely*fpair;



        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (newton_pair || j < nlocal) {

          local_force_x[atom_id[j]-1] -= delx*fpair;
          local_force_y[atom_id[j]-1] -= dely*fpair;
	      f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
         
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }
    MPI_Allreduce(local_force_x,force_x,n_atoms,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(local_force_y,force_y,n_atoms,MPI_DOUBLE,MPI_SUM,world);
    
    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

    //for (i = 1; i <= n_atoms; i++) {
        force_atom_x[i] = force_x[atom_id[i]-1];
        force_atom_y[i] = force_y[atom_id[i]-1];
        }
  if (vflag_fdotr)
    virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   This is used be pair_write;
   it is called only if rsq < cutsq[itype][jtype], no need to check that
------------------------------------------------------------------------- */

double PairWCAatt::single(int /* i */, int /* j */, int itype, int jtype, double rsq,
                         double /* factor_coul */, double factor_lj,
                         double &fforce)
{
  double r, r2inv, r6inv, cosone, force, energy;

  r = sqrt(rsq);
  force = energy = 0;
  if (r <= sigma[itype][jtype]) {
      r2inv = 1.0/rsq;
      r6inv = r2inv*r2inv*r2inv;
      force = r6inv*(lj12_f[itype][jtype]*r6inv - lj6_f[itype][jtype])*r2inv;
      energy = r6inv*(lj12_e[itype][jtype]*r6inv - lj6_e[itype][jtype]);
      energy += eps_wca[itype][jtype]-epsilon[itype][jtype];
  }else{
    cosone = cos(MY_PI*(r-sigma[itype][jtype]) / (2.0*w[itype][jtype]));
    force = -(MY_PI*epsilon[itype][jtype] / (2.0*w[itype][jtype])) *
                 sin(MY_PI*(r-sigma[itype][jtype]) / w[itype][jtype]) / r;
    energy = -epsilon[itype][jtype]*cosone*cosone;
  }
  fforce = factor_lj*force;
  return factor_lj*energy;
}


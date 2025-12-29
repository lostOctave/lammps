/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond_harmonic_area.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
// to check
#include "update.h"
#include <iostream>
using namespace LAMMPS_NS;
using namespace std;
/* ---------------------------------------------------------------------- */

AreaBondHarmonic::AreaBondHarmonic(LAMMPS *lmp) : Bond(lmp)
{
  reinitflag = 1;

  // chunk-based data
  nchunk = 1;
  maxchunk = 0;
  allocate();
}

/* ---------------------------------------------------------------------- */

AreaBondHarmonic::~AreaBondHarmonic()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(area);
    memory->destroy(area_all);
    memory->destroy(peri_local);
    memory->destroy(peri_all);
    memory->destroy(local_max_x);
    memory->destroy(local_min_x);
    memory->destroy(max_x);
    memory->destroy(min_x);
    memory->destroy(local_max_y);
    memory->destroy(local_min_y);
    memory->destroy(max_y);
    memory->destroy(min_y);
    memory->destroy(local_xc);
    memory->destroy(local_yc);
    memory->destroy(xc);
    memory->destroy(yc);
    memory->destroy(local_size_cell);
    memory->destroy(size_cell);
    memory->destroy(bond_left);
    memory->destroy(bond_right);
  }
}

/* ---------------------------------------------------------------------- */

void AreaBondHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond, sumx;
  double rsq,r,dr,rk;

  ebond = 0.0;
  ev_init(eflag,vflag);

//  if( update->ntimestep >= 50000){ printf("START COMPUTE\n"); }


  //Fix ->v_tally(int, double *);
  nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();
  int *ichunk = cchunk->ichunk;
  if (nchunk > maxchunk) allocate();

  for(int i = 0; i < nchunk; i++){
    area[i] = 0.0;
    peri_local[i] = 0.0;
    local_xc[i] = local_yc[i] = xc[i] = yc[i] = 0;
    local_size_cell[i] = size_cell[i] = 0;
    local_max_x[i] = domain->boxlo[0]-1;
    local_min_x[i] = domain->boxhi[0]+1;
    local_max_y[i] = domain->boxlo[1]-1;
    local_min_y[i] = domain->boxhi[1]+1;
  }

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int **bondlist = neighbor->bondlist;
  imageint *image = atom->image;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int n_atoms = atom->natoms;
  int newton_bond = force->newton_bond;
  int *atom_id = atom->tag;
  int **atom_special = atom->special;
// image and unwrap coordinates of the two atoms
  double ux_i1[3], ux_i2[3];

  int flag;
  int cols;

  int index1; //index2; // index of the chunks of the two atoms
  // we assume are equal

  char length_bond[] = "length_bond";
  int index_bond = atom->find_custom(length_bond,flag, cols);
  double *bond_lengths = atom->dvector[index_bond];

  char area_cell[] = "area_cell";
  int index_area = atom->find_custom(area_cell,flag, cols);
  double *a0 = atom->dvector[index_area];

  char karea_cell[] = "k_area";
  int index_karea = atom->find_custom(karea_cell,flag, cols);
  double *karea = atom->dvector[index_karea];

  char XC[] = "xc";
  int index_xc = atom->find_custom(XC,flag, cols);
  double *XCcell = atom->dvector[index_xc];

  char YC[] = "yc";
  int index_yc = atom->find_custom(YC,flag, cols);
  double *YCcell = atom->dvector[index_yc];

  char AR[] = "area";
  int index_ar = atom->find_custom(AR,flag,cols);
  double *area_store = atom->dvector[index_ar];

  char PERI[] = "peri";
  int index_peri = atom->find_custom(PERI,flag,cols);
  double *peri = atom->dvector[index_peri];

  char SI[] = "shape";
  int index_si = atom->find_custom(SI,flag,cols);
  double *q0 = atom->dvector[index_si];

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
  tagint *molecule = atom->molecule;/// Find molecule id
  // compute the minimum and max x of each particle;
  // this is needed to assess if a particle cross a bc when computing the area
  for(int i = 0; i < nlocal; i++){
    //index1 = ichunk[i]-1;
    index1 = molecule[i]-1;
    if(index1 >= 0){
      // extreme of the cell for area calculation
      if( x[i][0] > local_max_x[index1] ) local_max_x[index1] = x[i][0];
      if( x[i][0] < local_min_x[index1] ) local_min_x[index1] = x[i][0];
      if( x[i][1] > local_max_y[index1] ) local_max_y[index1] = x[i][1];
      if( x[i][1] < local_min_y[index1] ) local_min_y[index1] = x[i][1];
      // count the molecule in the chunk
      local_size_cell[index1]++;
      // unwrap coordinate for COM computation
      domain->unmap(x[i],image[i],unwrap);
      local_xc[index1] += unwrap[0];
      local_yc[index1] += unwrap[1];
    }
  }
  // extreme for area calculation
  MPI_Allreduce(local_max_x,max_x,nchunk,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(local_min_x,min_x,nchunk,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(local_max_y,max_y,nchunk,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(local_min_y,min_y,nchunk,MPI_DOUBLE,MPI_MIN,world);
  // center of mass
  MPI_Allreduce(local_xc,xc,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(local_yc,yc,nchunk,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(local_size_cell,size_cell,nchunk,MPI_INT,MPI_SUM,world);
  for(int i = 0; i < nchunk; i++){
    xc[i] /= size_cell[i];
    yc[i] /= size_cell[i];
  }

  for (n = 0; n < nbondlist; n++) {
    type = bondlist[n][2];
    if(type == 1){ // if this is an areabond
      i1 = bondlist[n][0];
      i2 = bondlist[n][1];
      domain->unmap(x[i1],image[i1],ux_i1);
      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];
      delz = x[i1][2] - x[i2][2];
      rsq = delx*delx + dely*dely; // + delz*delz;
      r = sqrt(rsq);
      
/////////////////////////////////////////////////////////////////////// each atom know the equilibrium length of its bond
//    dr = r - r0[type];
      dr = r - bond_lengths[i1];
      rk = k[type] * dr;
      // force & energy
      if (r > 0.0) fbond = -2.0*rk/r;
      else fbond = 0.0;
      if (eflag) ebond = rk*dr;

      // apply force to each of 2 atoms
      if (newton_bond || i1 < nlocal) {
        f[i1][0] += delx*fbond;
        f[i1][1] += dely*fbond;
        f[i1][2] += delz*fbond;
      }

      if (newton_bond || i2 < nlocal) {
        f[i2][0] -= delx*fbond;
        f[i2][1] -= dely*fbond;
        f[i2][2] -= delz*fbond;
      }

      // contribution of the bond to the area; we get the unwrapped coordinate of the first particle in the bond, and its chunk
      //index1 = ichunk[i1]-1;
      index1 = molecule[i1]-1;

      if(index1 >= 0){
        x1 = ux_i1[0];
        y1 = ux_i1[1];
        area[index1] += 0.5*(2*x1-delx)*(dely); /// (1/2)*sumx*dy

        peri_local[index1] += r;
      }
      if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
    }

  }

  MPI_Allreduce(area,area_all,nchunk,MPI_DOUBLE,MPI_SUM,world);

  MPI_Allreduce(peri_local,peri_all,nchunk,MPI_DOUBLE,MPI_SUM,world);

  // now, area_all is a vector containing the area of each chunk

//  error->all(FLERR,"check");
  double prefactor, f0, f1,v[6];
  //void v_tally = (Fix *) ->v_tally(int,double *);
  for (n = 0; n < nbondlist; n++) {
    type = bondlist[n][2];
    if(type == 1){
      i1 = bondlist[n][0];
      i2 = bondlist[n][1];

      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];


      domain->unmap(x[i1],image[i1],ux_i1);
      ux_i2[0] = ux_i1[0]-delx;
      ux_i2[1] = ux_i1[1]-dely;
      //index1 = ichunk[i1]-1;
      index1 = molecule[i1]-1;
      if(index1 >= 0){
	area_store[i1] = area_all[index1];

	peri[i1] = peri_all[index1];
	q0[i1] = peri[i1]/sqrt(area_store[i1]);
        prefactor = -karea[i1]*(area_all[index1]-a0[i1]);

        f0 = prefactor*(ux_i1[1]-ux_i2[1]); // 2*da/dx_i1 = 2*da/dx_i2
        f1 = prefactor*(ux_i2[0]+ux_i1[0]); // 2*da/dy_i1 = -2*da/dy_i2
        XCcell[i1] = xc[index1];
	YCcell[i1] = yc[index1];
	XCcell[i2] = xc[index1];
        YCcell[i2] = yc[index1];
	f[i1][0] += f0;
        f[i2][0] += f0;
        f[i1][1] += f1;
        f[i2][1] -= f1;
      }
     //if(eflag)  XCcell[i1] = eatom[i1];

        //YCcell[i1] = a0[i1];//molecule[i1];
      ebond = 0;

      if(eflag){ ebond = 0.5*karea[i1]*powl(area_all[index1]-a0[i1],2)/size_cell[index1]; }
      if(evflag){

        ev_area_tally(i1,ebond,f0,f1,ux_i1[0]-XCcell[i1],ux_i1[1]-YCcell[i1],ux_i1[2]);
        ev_area_tally(i2,ebond,f0,-f1,ux_i2[0]-XCcell[i2],ux_i2[1]-YCcell[i2],ux_i2[2]);
      }

    }
  }

}

void AreaBondHarmonic::ev_area_tally(int i,
                    double ebond, double f0, double f1,
                    double x, double y, double z)
{
  double ebondhalf,v[6];

          v[0] = f0 * x;
          v[1] = f1 * y;
          v[2] = 0 * z;
          v[3] = f0 * y;
          v[4] = f0 * z;
          v[5] = f1 * z;
  if (eflag_either) {
    if (eflag_global) {
      energy += ebond;
    }
    if (eflag_atom) {
      eatom[i] += ebond;
    }
  }
if (vflag_either) {

    if (vflag_global) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
     }
   if (vflag_atom) {
	vatom[i][0] += v[0];
        vatom[i][1] += v[1];
        vatom[i][2] += v[2];
        vatom[i][3] += v[3];
        vatom[i][4] += v[4];
        vatom[i][5] += v[5];
	}
    }
}


/* ---------------------------------------------------------------------- */
void AreaBondHarmonic::allocate_bonds()
{
    allocated = 1;
    int n = atom->nbondtypes;
    memory->create(k,n+1,"bond:k");
    memory->create(r0,n+1,"bond:r0");
    memory->create(setflag,n+1,"bond:setflag");
    for (int i = 1; i <= n; i++) setflag[i] = 0;
}

void AreaBondHarmonic::allocate()
{
  // compute chunk/atom assigns atoms to chunk IDs
  // extract ichunk index vector from compute
  // ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

  // allocate the chunks
  maxchunk = nchunk;
  //int n_atoms = atom->natoms;
  int n_local = atom->nlocal + atom->nghost;
//  printf("Num chunks: %d\n",nchunk);
  memory->create(area,maxchunk,"bond:area");
  memory->create(area_all,maxchunk,"bond:area_all");
  memory->create(peri_all,maxchunk,"bond:peri_all");
  
  
  memory->create(peri_local,maxchunk,"bond:peri_local");
  memory->create(local_max_x,maxchunk,"bond:local_max_x");
  memory->create(local_min_x,maxchunk,"bond:local_min_x");
  memory->create(max_x,maxchunk,"bond:max_x");
  memory->create(min_x,maxchunk,"bond:min_x");
  memory->create(local_max_y,maxchunk,"bond:local_max_y");
  memory->create(local_min_y,maxchunk,"bond:local_min_y");
  memory->create(max_y,maxchunk,"bond:max_y");
  memory->create(min_y,maxchunk,"bond:min_y");

  memory->create(local_xc,maxchunk,"bond:local_xc");
  memory->create(xc,maxchunk,"bond:xc");
  memory->create(local_yc,maxchunk,"bond:local_yc");
  memory->create(yc,maxchunk,"bond:yc");

  memory->create(local_size_cell,maxchunk,"bond:local_size_cell");
  memory->create(size_cell,maxchunk,"bond:size_cell");
  memory->create(bond_left,n_local+1,2,"bond:bond_left");
  memory->create(bond_right,n_local+1,2,"bond:bond_right");
}



/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AreaBondHarmonic::coeff(int narg, char **arg)
{

  if (narg != 3) error->all(FLERR,"Incorrect args for bond coefficients: k r0 chunk_id");
  if (!allocated) allocate_bonds();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double k_one = utils::numeric(FLERR,arg[1],false,lmp);
  double r0_one = utils::numeric(FLERR,arg[2],false,lmp);

  sprintf(arg[0],"Cells");
  int n = strlen(arg[0]) + 1;
  idchunk = new char[n];
  strcpy(idchunk,arg[0]);
  int icompute = modify->find_compute(idchunk);
  if (icompute >= 0) {
    cchunk = (ComputeChunkAtom *) modify->compute[icompute];
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double AreaBondHarmonic::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AreaBondHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AreaBondHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&r0[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AreaBondHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],r0[i]);
}

/* ---------------------------------------------------------------------- */

double AreaBondHarmonic::single(int type, double rsq, int /*i*/, int /*j*/,
                            double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr;
  fforce = 0;
  if (r > 0.0) fforce = -2.0*rk/r;
  return rk*dr;
}

/* ----------------------------------------------------------------------
    Return ptr to internal members upon request.
------------------------------------------------------------------------ */
void *AreaBondHarmonic::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str,"kappa")==0) return (void*) k;
  if (strcmp(str,"r0")==0) return (void*) r0;
  return nullptr;
}

double AreaBondHarmonic::findCircle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double x12 = x1 - x2;
    double x13 = x1 - x3;

    double y12 = y1 - y2;
    double y13 = y1 - y3;

    double y31 = y3 - y1;
    double y21 = y2 - y1;

    double x31 = x3 - x1;
    double x21 = x2 - x1;

    // x1^2 - x3^2
    double sx13 = pow(x1, 2) - pow(x3, 2);

    // y1^2 - y3^2
    double sy13 = pow(y1, 2) - pow(y3, 2);

    double sx21 = pow(x2, 2) - pow(x1, 2);
    double sy21 = pow(y2, 2) - pow(y1, 2);

    double f = ((sx13) * (x12)
             + (sy13) * (x12)
             + (sx21) * (x13)
             + (sy21) * (x13))
            / (2 * ((y31) * (x12) - (y21) * (x13)));
    double g = ((sx13) * (y12)
             + (sy13) * (y12)
             + (sx21) * (y13)
             + (sy21) * (y13))
            / (2 * ((x31) * (y12) - (x21) * (y13)));

    double c = -pow(x1, 2) - pow(y1, 2) - 2 * g * x1 - 2 * f * y1;

    // eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
    // where centre is (h = -g, k = -f) and radius r
    // as r^2 = h^2 + k^2 - c
    double h = -g;
    double k = -f;
    double sqr_of_r = h * h + k * k - c;

    // r is the radius
    double r = sqrt(sqr_of_r);

    //cout << "Centre = (" << h << ", " << k << ")" << endl;
    //cout << "Radius = " << r;
    return r;
}

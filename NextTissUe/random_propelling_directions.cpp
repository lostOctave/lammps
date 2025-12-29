#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "xrand.h"
#include "iostream"
#include "fstream"

using namespace std;


int main(int narg, char **arg){
  int NL,NS,nbL,nbS;
  sscanf(arg[1],"%d",&NL);
  sscanf(arg[2],"%d",&nbL);
  sscanf(arg[3],"%d",&NS);
  sscanf(arg[4],"%d",&nbS);

  
  int id = 0;
  int ntot =  NL*nbL+NS*nbS;
  char buffer[100];
  sprintf(buffer,"tail -n %d mol_id.txt > idmol.txt",ntot);
  system(buffer);
  FILE* in = fopen("idmol.txt","r");
  int mol;
  ofstream outfile;
  outfile.open("self_propel_theta.txt");
  int nmol = NL+NS;
  double theta[nmol];
  for(int i = 0; i < NL+NS; i++) theta[i] =  Xrandom()*2*M_PI;
  outfile << ntot << endl;
  while(fscanf(in,"%d %d",&id, &mol) > 0){
    if(mol != 0)
      outfile << id << " " << theta[mol-1] << endl;
  }
  outfile.close();
}

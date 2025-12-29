#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "iostream"
#include "fstream"

using namespace std;


int main(int narg, char **arg){
  int N;
  double Radius;
  sscanf(arg[1],"%d",&N);
  sscanf(arg[2],"%lf",&Radius);
  double circ = N;
  char line[1024];
//  system("rm template.cell"); 
  sprintf(line,"template.cell");
  ofstream outfile(line);
  outfile << "#template for a cell with " << N  << " atoms " << endl << endl;
  outfile << N << " atoms" << endl;
  outfile << N << " bonds" << endl;
//  outfile << N-2 << " angles" << endl;
  outfile << endl;
  outfile << "Coords" << endl << endl;
  for(int n = 0; n < N; n++){
    outfile << n+1 << "	" << Radius*cos(n*2*M_PI/N) << "	" << -Radius*sin(n*2*M_PI/N) << "	" << 0 <<  endl; // id x y z
  }
  // central particle 
//  outfile << N+1 << "	" << 0 << "	" << 0 << "	" << 0 << endl;
  outfile << endl;
  outfile << "Molecules" << endl << endl;
  for(int n = 0; n < N; n++){
    outfile << n+1 << "	" << 1 << endl; // id mol-id
  }
  // central particle 
//  outfile << N+1 << "	" << 1 << "	" << endl;
  outfile << endl;

  outfile << "Types" << endl << endl;
  for(int n = 0; n < N; n++) outfile << n+1 << "	" << 1 << endl; // border particles are of type 1
//  outfile << N+1 << "	" << 2 << endl; // central particle is of type 2

  outfile << endl;
  outfile << "Bonds" << endl << endl;
  for(int n = 1; n <= N; n++){
    if(n < N)
      outfile << n << "	" << 1 << "	" << n << "	" << n+1 << endl;
    else
      outfile << n << "	" << 1 << "	" << n << "	" << 1 << endl;
  }
  outfile << endl;
//  outfile << "Angles" << endl << endl;
//  for(int n = 2; n < N; n++) outfile << n-1 << "	" << 1 << "	" << n-1 << "	" << n << "	" << n+1 << endl;
  outfile.close();
}


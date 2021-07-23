#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char **argv){
  unsigned int myid, numprocs, icp = 4, n = 100;
  unsigned int s, ts, ia, na, *a;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  na = n / icp;
  if((a = (int *)calloc(na, sizeof(int))) == NULL){
    fprintf(stderr, "malloc error\n");
    return 1;
  }

  int ist = myid * na;
  int iet = (myid + 1) * na;
  for(int i = ist; i < iet; i++){
    ia = i - (myid * na);
    a[ia] = i + 1;
  }
  if(myid == 1){
    for(int i = ist; i < iet; i++) fprintf(stderr, "a[%d] = %d\n", i, a[i]);
  }

  s = 0;
  ts = 0;
  for(int i = 0; i < na; i++) s += a[i];

  MPI_Reduce(&s, &ts, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  fprintf(stderr, "irank = %d, s = %d, ts = %d\n", myid, s, ts);

  MPI_Finalize();

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char **argv){
  int *a, n, s, myid, numprocs, ierr;
  int ts, ipe, ist, iet;

  //fprintf(stderr, "Input num : ");
  //scanf("%d", &n);
  n = 100;
  if((a = (int *)calloc(n, sizeof(int))) == NULL){
    fprintf(stderr, "allocation error\n");
    return 1;
  }

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  fprintf(stderr, "myid = %d, numprocs = %d\n", myid, numprocs);

  for(int i = 0; i < n; i++) a[i] = i;

  s = 0;
  ts = 0;
  ipe = n / numprocs;
  ist = 1 + myid * ipe;
  iet = (myid + 1) * ipe;

  for(int i = ist; ist < iet; ist++) s += a[i];

  MPI_Reduce(&s, &ts, 1, MPI_INT, MPI_SUM, myid+1, MPI_COMM_WORLD);
  fprintf(stderr, "myid = %d %d %d\n", myid, s, ts);

  ierr = MPI_Finalize();

  return 0;
}

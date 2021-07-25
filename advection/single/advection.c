#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allocate.h"

void input(int*, int*, double*, double*, double*, double*, double*, int*, int*, double*, double*, int**, int*, int*, double*, double*, double*, double*, double***);

void init(int*, double*, double*, double*, int*, int*, double*);

void matrix(int*, int**, double*, double*, double*, double*, double*, double*, double*, double***, double*);

void vec(int*, int*, int**, double***, double*, double*);

void solve(int*, double*, double*);

void bound(int*, int*, double*, double*);

void output(int*, double*, int, double*, int*, int**, double*, double*);

void change(int*, double*, double*);

int main(int argc, char **argv){
  double *xx, *yy, *fbc, ***rmat, *amb, *cc1, *cc2;
  int **nc, *nbc;

  int itmax, iout, nnod, nelem, ibc;
  double uu, vv, skx, sky, dt;

  FILE *fp1, *fp2, *fp3;
  int j;

  fp1 = fopen("test/input.dat", "r");

  if(fp1 == NULL){
    fprintf(stderr, "Can't open file\n");
    exit(1);
  }

  //fprintf(stderr, "flag1\n");
  fscanf(fp1, "%d %d %lf", &itmax, &iout, &dt);
  //fprintf(stderr, "flag2\n");
  fscanf(fp1, "%lf %lf %lf %lf", &uu, &vv, &skx, &sky);
  //fprintf(stderr, "%d %d %lf\n", *itmax, *iout, *dt);
  //fprintf(stderr, "%lf %lf %lf %lf\n", *uu, *vv, *skx, *sky);

  fclose(fp1);

  fp2 = fopen("test/mesh.dat", "r");

  //fprintf(stderr, "flag2\n");
  fscanf(fp2, "%d %d", &nnod, &nelem);

  xx = (double*)allocate_vector(sizeof(double), nnod);
  yy = (double*)allocate_vector(sizeof(double), nnod);
  nc = (int**)allocate_matrix(sizeof(int), 3, nelem);
  amb = (double*)allocate_vector(sizeof(double), nnod);
  cc1 = (double*)allocate_vector(sizeof(double), nnod);
  cc2 = (double*)allocate_vector(sizeof(double), nnod);
  rmat = (double***)malloc(3 * sizeof(double));
  for(int i = 0; i < 3; i++){
    rmat[i] = (double**)malloc(3 * sizeof(double));
    for(int j = 0; j < 3; j++) rmat[i][j] = (double*)malloc(nelem * sizeof(double));
  }

  for(int i = 0; i < nnod; i++) fscanf(fp2, "%d %lf %lf", &j, &xx[i], &yy[i]);
  //for(int i = 0; i < *nnod; i++) fprintf(stderr, "%d %lf %lf\n", i, xx[i], yy[i]);
  for(int i = 0; i < nelem; i++) fscanf(fp2, "%d %d %d %d", &j, &nc[0][i], &nc[1][i], &nc[2][i]);
  for(int i = 0; i < nelem; i++){
    nc[0][i] -= 1;
    nc[1][i] -= 1;
    nc[2][i] -= 1;
    //fprintf(stderr, "%d %d %d %d\n", i, nc[0][i], nc[1][i], nc[2][i]);
  }

  fclose(fp2);

  fp3 = fopen("test/bc.dat", "r");

  //fprintf(stderr, "flag3\n");
  fscanf(fp3, "%d", &ibc);
  if(ibc != 0){    
    nbc = (int*)allocate_vector(sizeof(int), ibc);
    fbc = (double*)allocate_vector(sizeof(double), ibc);
    for(int i = 0; i < ibc; i++) fscanf(fp3, "%d %d %lf", &j, &nbc[i], &fbc[i]);
    for(int i = 0; i < ibc; i++) nbc[i] -= 1;
    //for(int i = 0; i < *ibc; i++) fprintf(stderr, "nbc[%d] = %d\n", i, nbc[i]);
  }

  fclose(fp3);

  for(int i = 0; i < nnod; i++){
    //fprintf(stderr, "%d\n", i);
    cc1[i] = 0.0;
    //fprintf(stderr, "%d\n", i);
    cc2[i] = 0.0;
    //fprintf(stderr, "%d\n", i);
    amb[i] = 0.0;
    //fprintf(stderr, "%d\n", i);
    //fprintf(stderr, "cc1[%d] = %lf cc2[%d] = %lf amb[%d] = %lf\n", i, cc1[i], i, cc2[i], i, amb[i]);
  }

  int ii;
  for(int i = 0; i < ibc; i++){
    ii = nbc[i];
    cc1[ii] = fbc[i];
  }

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      for(int im = 0; im < nelem; im++) rmat[i][j][im] = 0.0;
    }
  }

  fprintf(stderr, "input\n");
  //input(&itmax, &iout, &dt, &uu, &vv, &skx, &sky, &nnod, &nelem, xx, yy, nc, &ibc, nbc, fbc, amb, cc1, cc2, rmat);

  fprintf(stderr, "init\n");
  //init(&nnod, cc1, cc2, amb, &ibc, nbc, fbc);

  fprintf(stderr, "matrix\n");
  matrix(&nelem, nc, xx, yy, amb, &uu, &vv, &skx, &sky, rmat, &dt);

  for(int istep = 0; istep <= itmax; istep++){
    fprintf(stderr, "vec\n");
    vec(&nnod, &nelem, nc, rmat, cc1, cc2);

    fprintf(stderr, "solve\n");
    solve(&nnod, amb, cc2);

    fprintf(stderr, "bound\n");
    bound(&ibc, nbc, fbc, cc2);

    if(istep % iout == 0) output(&nnod, cc2, istep, &dt, &nelem, nc, xx, yy);

    fprintf(stderr, "change\n");
    change(&nnod, cc1, cc2);
  }

  fprintf(stderr, "xx\n");
  free(xx);
  fprintf(stderr, "yy\n");
  free(yy);
  fprintf(stderr, "fbc\n");
  free(fbc);
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++) free(rmat[i][j]);
    free(rmat[i]);
  }
  fprintf(stderr, "rmat\n");
  free(rmat);
  fprintf(stderr, "amb\n");
  free(amb);
  fprintf(stderr, "cc1\n");
  free(cc1);
  fprintf(stderr, "cc2\n");
  free(cc2);
  for(int i = 0; i < 3; i++) free(nc[i]);
  fprintf(stderr, "nc\n");
  free(nc);
  fprintf(stderr, "nbc\n");
  free(nbc);

  return 0;
}


void input(int* itmax, int* iout, double* dt, double* uu, double* vv, double* skx, double* sky, int* nnod, int* nelem, double* xx, double* yy, int** nc, int* ibc, int* nbc, double* fbc, double* amb, double* cc1, double* cc2, double ***rmat){
  FILE *fp1, *fp2, *fp3;
  int j;

  fp1 = fopen("test/input.dat", "r");

  if(fp1 == NULL){
    fprintf(stderr, "Can't open file\n");
    exit(1);
  }

  //fprintf(stderr, "flag1\n");
  fscanf(fp1, "%d %d %lf", itmax, iout, dt);
  //fprintf(stderr, "flag2\n");
  fscanf(fp1, "%lf %lf %lf %lf", uu, vv, skx, sky);
  //fprintf(stderr, "%d %d %lf\n", *itmax, *iout, *dt);
  //fprintf(stderr, "%lf %lf %lf %lf\n", *uu, *vv, *skx, *sky);

  fclose(fp1);

  fp2 = fopen("test/mesh.dat", "r");

  //fprintf(stderr, "flag2\n");
  fscanf(fp2, "%d %d", nnod, nelem);

  xx = (double*)allocate_vector(sizeof(double), *nnod);
  yy = (double*)allocate_vector(sizeof(double), *nnod);
  nc = (int**)allocate_matrix(sizeof(int), 3, *nelem);
  amb = (double*)allocate_vector(sizeof(double), *nnod);
  cc1 = (double*)allocate_vector(sizeof(double), *nnod);
  cc2 = (double*)allocate_vector(sizeof(double), *nnod);
  rmat = (double***)malloc(3 * sizeof(double));
  for(int i = 0; i < 3; i++){
    rmat[i] = (double**)malloc(3 * sizeof(double));
    for(int j = 0; j < 3; j++) rmat[i][j] = (double*)malloc(*nelem * sizeof(double));
  }

  for(int i = 0; i < *nnod; i++) fscanf(fp2, "%d %lf %lf", &j, &xx[i], &yy[i]);
  //for(int i = 0; i < *nnod; i++) fprintf(stderr, "%d %lf %lf\n", i, xx[i], yy[i]);
  for(int i = 0; i < *nelem; i++) fscanf(fp2, "%d %d %d %d", &j, &nc[0][i], &nc[1][i], &nc[2][i]);
  for(int i = 0; i < *nelem; i++){
    nc[0][i] -= 1;
    nc[1][i] -= 1;
    nc[2][i] -= 1;
    //fprintf(stderr, "%d %d %d %d\n", i, nc[0][i], nc[1][i], nc[2][i]);
  }

  fclose(fp2);

  fp3 = fopen("test/bc.dat", "r");

  //fprintf(stderr, "flag3\n");
  fscanf(fp3, "%d", ibc);
  if(*ibc != 0){    
    nbc = (int*)allocate_vector(sizeof(int), *ibc);
    fbc = (double*)allocate_vector(sizeof(double), *ibc);
    for(int i = 0; i < *ibc; i++) fscanf(fp3, "%d %d %lf", &j, &nbc[i], &fbc[i]);
    for(int i = 0; i < *ibc; i++) nbc[i] -= 1;
    //for(int i = 0; i < *ibc; i++) fprintf(stderr, "nbc[%d] = %d\n", i, nbc[i]);
  }

  fclose(fp3);

  for(int i = 0; i < *nnod; i++){
    //fprintf(stderr, "%d\n", i);
    cc1[i] = 0.0;
    //fprintf(stderr, "%d\n", i);
    cc2[i] = 0.0;
    //fprintf(stderr, "%d\n", i);
    amb[i] = 0.0;
    //fprintf(stderr, "%d\n", i);
    fprintf(stderr, "cc1[%d] = %lf cc2[%d] = %lf amb[%d] = %lf\n", i, cc1[i], i, cc2[i], i, amb[i]);
  }

  int ii;
  for(int i = 0; i < *ibc; i++){
    ii = nbc[i];
    cc1[ii] = fbc[i];
  }

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      for(int im = 0; im < *nelem; im++) rmat[i][j][im] = 0.0;
    }
  }
}

void init(int* nnod, double* cc1, double* cc2, double* amb, int* ibc, int* nbc, double* fbc){
  int ii;

  fprintf(stderr, "flag1\n");
  for(int i = 0; i < *nnod; i++){
    fprintf(stderr, "%d\n", i);
    cc1[i] = 0.0;
    fprintf(stderr, "%d\n", i);
    cc2[i] = 0.0;
    fprintf(stderr, "%d\n", i);
    amb[i] = 0.0;
    fprintf(stderr, "%d\n", i);
  }
  fprintf(stderr, "flag2\n");
  for(int i = 0; i < *ibc; i++){
    ii = nbc[i];
    cc1[ii] = fbc[i];
  }
  fprintf(stderr, "flag3\n");
}

void matrix(int* nelem, int** nc, double* xx, double* yy, double* amb, double* uu, double* vv, double* skx, double* sky, double*** rmat, double* dt){
  int ia, ib, ic;
  double xa, xb, xc, ya, yb, yc;
  double area, area1, area2, *b, *c, aa06, aa12, aa03;

  b = (double*)allocate_vector(sizeof(double), 3);
  c = (double*)allocate_vector(sizeof(double), 3);

  /*for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      for(int im = 0; im < *nelem; im++) rmat[i][j][im] = 0.0;
    }
  }*/

  //fprintf(stderr, "flag1\n");
  for(int im = 0; im < *nelem; im++){
    ia = nc[0][im];
    ib = nc[1][im];
    ic = nc[2][im];
    //fprintf(stderr, "flag2\n");
    xa = xx[ia];
    xb = xx[ib];
    xc = xx[ic];
    ya = yy[ia];
    yb = yy[ib];
    yc = yy[ic];
    //fprintf(stderr, "flag3\n");
    area = (xa * (yb - yc) + xb * (yc - ya) + xc * (ya - yb)) * 0.5;
    area1 = area * 2.0;
    area2 = 1.0 / area1;
    //fprintf(stderr, "flag4\n");
    if(area <= 0.0){
      fprintf(stderr, "area <= 0.0\n");
      exit(1);
    }
    //fprintf(stderr, "flag5\n");
    b[0] = (yb - yc) * area2;
    b[1] = (yc - ya) * area2;
    b[2] = (ya - yb) * area2;
    c[0] = (yc - yb) * area2;
    c[1] = (ya - yc) * area2;
    c[2] = (yb - ya) * area2;

    aa06 = area / 6.0;
    aa12 = area / 12.0;

    rmat[0][0][im] = aa06;
    rmat[0][1][im] = aa12;
    rmat[0][2][im] = aa12;
    rmat[1][0][im] = aa12;
    rmat[1][1][im] = aa06;
    rmat[1][2][im] = aa12;
    rmat[2][0][im] = aa12;
    rmat[2][1][im] = aa12;
    rmat[2][2][im] = aa06;

    aa03 = area / 3.0;
    //fprintf(stderr, "flag6\n");
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < 3; j++){
        rmat[i][j][im] -= (*dt) * ((*uu) * aa03 * b[j] + (*vv) * aa03 * c[j]);
      }
    }
    //fprintf(stderr, "flag7\n");
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < 3; j++){ 
        rmat[i][j][im] -= (*dt) * (b[i] * b[j] * (*skx) * area + c[i] * c[j] * (*sky) * area);
      } 
    }
    //fprintf(stderr, "flag8\n");
    amb[ia] += aa03;
    amb[ib] += aa03;
    amb[ic] += aa03;   
  }

  free(b);
  free(c);
}

void vec(int* nnod, int* nelem, int** nc, double*** rmat, double* cc1, double* cc2){
  int ii, jj;

  for(int i = 0; i < *nnod; i++) cc2[i] = 0.0;

  for(int im = 0; im < *nelem; im++){
    for(int i = 0; i < 3; i++){
      ii = nc[i][im];
      for(int j = 0; j < 3; j++){
        jj = nc[j][im];
	cc2[ii] += rmat[i][j][im] * cc1[jj];
      }
    }
  }
}

void solve(int* nnod, double* amb, double* cc2){
  double ss;

  for(int i = 0; i < *nnod; i++){
    ss = 1.0 / amb[i];
    cc2[i] = ss * cc2[i];
  }
}

void bound(int* ibc, int* nbc, double* fbc, double* cc2){
  int ii;

  for(int i = 0; i < *ibc; i++){
    ii = nbc[i];
    cc2[ii] = fbc[i];
  }
}

void output(int* nnod, double* cc2, int istep, double* dt, int* nelem, int** nc, double* xx, double* yy){
  FILE *fp;
  char filepath[256];
  double time = (double)istep * (*dt);

  /*fprintf(stderr, "time = %lf\n", time);
  for(int i = 0; i < *nnod; i++){
    fprintf(stderr, "%d %lf", i, cc2[i]);
  }*/

  sprintf(filepath, "result/advection%d.vtk", istep);

  fp = fopen(filepath, "w");

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "triangle\nASCII\nDATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS  %d  float\n", *nnod);
  for(int i = 0; i < *nnod; i++) fprintf(fp, "%lf %lf 0.0\n", xx[i], yy[i]);
  fprintf(fp, "CELLS %d  %d\n", *nelem, *nelem * 4);
  for(int i = 0; i < *nelem; i++) fprintf(fp, "3  %d  %d  %d\n", nc[0][i], nc[1][i], nc[2][i]);
  fprintf(fp, "CELL_TYPES %d\n", *nelem);
  for(int i = 0; i < *nelem; i++) fprintf(fp, "5\n");
  fprintf(fp, "POINT_DATA %d\n", *nnod);
  fprintf(fp, "SCALARS  cc2  float\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for(int i = 0; i < *nnod; i++) fprintf(fp, "%lf\n", cc2[i]);

  fclose(fp);

}

void change(int* nnod, double* cc1, double* cc2){
  for(int i = 0; i < *nnod; i++){
    cc1[i] = cc2[i];
    cc2[i] = 0.0;
  }
}
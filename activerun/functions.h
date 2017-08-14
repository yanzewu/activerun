#include <math.h>
#include <stdio.h>
#include <fstream>


/*********************FUNCTIONS**************************/

#include <new>

template<class atype> atype *create_1d_array(int);
template<class atype> atype *create_2d_array(int);

/******************FUNCTIONS FOR 1D ARRAY*****************/

template<class atype> 
atype *create_1d_array(int m1) 
{
  atype *array;
  try {
    array = new atype[m1];
  } catch (...) {
    printf("something uknown wrong with 1D array creation\n");
  }
  return array;
}

/**************FUNCTIONS FOR 2D ARRAY*********************/

template<class atype> 
atype **create_2d_array(int m1, int m2) 
{
  atype **array;
  try {
    array = new atype*[m1];
    for (int i=0; i<m1; i++) array[i] = new atype[m2];
  }  catch (...) {
    printf("something uknown wrong with 2D array creation\n");
  }
  return array;
}


/**************FUNCTIONS FOR 3D ARRAY*********************/

template<class atype> 
atype ***create_3d_array(int m1, int m2, int m3) 
{
  atype ***array;
  try {
    array = new atype**[m1];
    for (int i=0; i<m1; i++) {
      array[i] = new atype*[m2];
      for (int k=0; k<m2; k++) {
	array[i][k] = new atype[m3];
      }
    }
  }  catch (...) {
    printf("something uknown wrong with 3D array creation\n");
  }
  return array;
}

/*******************RANDOM # GENERATOR*********************/

double fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);

}

/******************LAMMPS RANDOM # GENERATOR****************/
//Marsaglia random number genertor

double unimars(int seed)
{
  int ij,kl,i,j,k,l,ii,jj,m, i97, j97;
  double s,t, c, cd, cm, *u;
  
  u = new double[97+1];
  
  ij = (seed-1)/30082;
  kl = (seed-1) - 30082*ij;
  i = (ij/177) % 177 + 2;
  j = ij %177 + 2;
  k = (kl/169) % 178 + 1;
  l = kl % 169;
  for (ii = 1; ii <= 97; ii++) {
    s = 0.0;
    t = 0.5;
    for (jj = 1; jj <= 24; jj++) {
      m = ((i*j) % 179)*k % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;
      if ((l*m) % 64 >= 32) s = s + t;
      t = 0.5*t;
    }
    u[ii] = s;
  }
  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;
  i97 = 97;
  j97 = 33;

  double uni = u[i97] - u[j97];
  if (uni < 0.0) uni += 1.0;
  u[i97] = uni;
  i97--;
  if (i97 == 0) i97 = 97;
  j97--;
  if (j97 == 0) j97 = 97;
  c -= cd;
  if (c < 0.0) c += cm;
  uni -= c;
  if (uni < 0.0) uni += 1.0;
  return uni;

}

/********************POINT IN TRIANGLE**********************/
//Barycentric coordinate system
int intersect(double x1, double y1, double x2, double y2, double x3, double y3, double xp, double yp)
{  
  int inside = 0;
  double detT = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);
  double L1 = ((y2-y3)*(xp-x3) + (x3-x2)*(yp-y3))/detT;
  double L2 = ((y3-y1)*(xp-x3) + (x1-x3)*(yp-y3))/detT;
  double L3 = 1-L1-L2;
  if (0<L1 && L1<1 && 0<L2 && L2<1 && 0<L3 && L3<1)
    inside = 1;
  return inside;
}

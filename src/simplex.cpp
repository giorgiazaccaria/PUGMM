/*
 * Program: simplex.cpp
 * Based on:
 * Author : Michael F. Hutt
 * http://www.mikehutt.com
 * 11/3/97
 *
 * Copyright (c) 1997-2004 <Michael F. Hutt>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 * An implementation of the Nelder-Mead simplex method.
 */

/*
 * Modified by Alexa A. Sochaniwsky on 2025-08-04 for integration with Rcpp.
 */

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>

#define MAX_IT      1000      /* maximum number of iterations */
#define ALPHA       1.0       /* reflection coefficient */
#define BETA        0.5       /* contraction coefficient */
#define GAMMA       2.0       /* expansion coefficient */

using namespace Rcpp;

// Function to convert R matrix to 2D array
double** matrix_to_carray(NumericMatrix X) {
  int nrow = X.nrow();
  int ncol = X.ncol();

  double** array = new double*[nrow];
  for(int i = 0; i < nrow; i++) {
    array[i] = new double[ncol];
    for(int j = 0; j < ncol; j++) {
      array[i][j] = X(i,j);
    }
  }
  return array;
}

// Function to free 2D array
void free_carray(double** array, int nrow) {
  for(int i = 0; i < nrow; i++) {
    delete[] array[i];
  }
  delete[] array;
}

// Wrapper function to call R function from C++
  double call_r_function(Function func, int n1, int p, double* v, double* index, double** X, double* z_g, double* mu_g, double** sigma) {
    NumericVector index_vec(index, index + p);
    NumericVector z_vec(z_g, z_g + n1);
    NumericVector mu_vec(mu_g, mu_g + p);
    NumericVector v_vec(v, v + p);
    NumericMatrix X_mat = NumericMatrix(n1, p);
    NumericMatrix sigma_mat = NumericMatrix(p, p);
    for(int i = 0; i < n1; i++) {
      for(int j = 0; j < p; j++) {
        X_mat(i,j) = X[i][j];
      }
    }

    for(int i = 0; i < p; i++) {
      for(int j = 0; j < p; j++) {
        sigma_mat(i,j) = sigma[i][j];
      }
    }

    NumericVector result = func(v_vec, index_vec, X_mat, z_vec, mu_vec, sigma_mat);
    return result[0];
  }

// [[Rcpp::export]]
NumericVector simplex(Function r_func,
                      int n1,
                      int p,
                      NumericVector index_vec,
                      NumericMatrix X_mat,
                      NumericVector z_vec,
                      NumericVector mu_vec,
                      NumericMatrix sigma_mat,
                      NumericVector start_vec,
                      double EPSILON,
                      double scale) {

  // Convert R objects to C-style
  double* index = index_vec.begin();
  double* z_g = z_vec.begin();
  double* mu_g = mu_vec.begin();
  double* start = start_vec.begin();
  double** X = matrix_to_carray(X_mat);
  double** sigma = matrix_to_carray(sigma_mat);

  // Original variables
  int vs;
  int vh;
  int vg;

  int i, j, m, row, n;
  //int k;
  int itr;


  double** v = new double*[p+1];
  for(i = 0; i <= p; i++){
    v[i] = new double[p];
  }

  double* f = new double[p+1];
  double* vr = new double[p];
  double* ve = new double[p];
  double* vc = new double[p];
  double* vm = new double[p];

  double fr, fe, fc;//, min;
  double fsum, favg, s, cent;

  n = 0;
  for(j=0;j<p;j++){
    if(index[j]==1){
      n +=1;
    }
  }

  double pn = scale*(sqrt(n+1)-1+n)/(n*sqrt(2));
  double qn = scale*(sqrt(n+1)-1)/(n*sqrt(2));

  for (i=0;i<n;i++) {
    v[0][i] = start[i];
  }

  for (i=1;i<=n;i++) {
    for (j=0;j<n;j++) {
      if (i-1 == j) {
	      v[i][j] = pn + start[j];
      }
      else {
	      v[i][j] = qn + start[j];
      }
    }
  }

  for (j = 0; j <= n; j++){
    f[j] = call_r_function(r_func, n1, p, v[j], index, X, z_g, mu_g, sigma);
  }

  //k = n + 1;

   /* begin the main loop of the minimization */
  for (itr = 1; itr <= MAX_IT; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vg]) {
	      vg = j;
      }
    }

    /* find the index of the smallest value */
    vs=0;
    for (j=0;j<=n;j++) {
      if (f[j] < f[vs]) {
	      vs = j;
      }
    }

    /* find the index of the second largest value */
    vh=vs;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
	      vh = j;
      }
    }

    /* calculate the centroid */
    for (j=0;j<=n-1;j++) {
      cent=0.0;
      for (m=0;m<=n;m++) {
	      if (m!=vg) {
	        cent += v[m][j];
	      }
      }
      vm[j] = cent/n;
    }

    /* reflect vg to new vertex vr */
    for (j=0;j<=n-1;j++) {
      vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];
    }
    fr = call_r_function(r_func, n1, p, vr, index, X, z_g, mu_g, sigma);
    //k++;

    /* added <= */
    if (fr <= f[vh] && fr > f[vs]) {
      for (j=0;j<=n-1;j++) {
	      v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }

    /* investigate a step further in this direction */
    /* added <= */
    if ( fr <=  f[vs]) {
      for (j=0;j<=n-1;j++) {
	      ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];
      }
      fe = call_r_function(r_func, n1, p, ve, index, X, z_g, mu_g, sigma);
      //k++;

      if (fe < fr) {
        for (j = 0; j < n; j++){
          v[vg][j] = ve[j];
        }
        f[vg] = fe;
      } else {
        for (j=0;j<=n-1;j++) {
	        v[vg][j] = vr[j];
	      }
	      f[vg] = fr;
      }
    }

    /* check to see if a contraction is necessary */
    if (fr > f[vh]) {
      for (j=0;j<=n-1;j++) {
	      vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];
      }
      fc = call_r_function(r_func, n1, p, vc, index, X, z_g, mu_g, sigma);
      //k++;

      if (fc < f[vg]) {
	      for (j=0;j<=n-1;j++) {
	        v[vg][j] = vc[j];
	      }
	      f[vg] = fc;
      }
      /* at this point the contraction is not successful,
	    we must halve the distance from vs to all the
	    vertices of the simplex and then continue.
	    10/31/97 - modified to account for ALL vertices.
      */
      else {
	      for (row=0;row<=n;row++) {
	        if (row != vs) {
	          for (j=0;j<=n-1;j++) {
	            v[row][j] = v[vs][j]+(v[row][j]-v[vs][j])/2.0;
	          }
	        }
	      }
        f[vg] = call_r_function(r_func, n1, p, v[vg], index, X, z_g, mu_g, sigma);
        //k++;
        f[vh] = call_r_function(r_func, n1, p, v[vh], index, X, z_g, mu_g, sigma);
        //k++;
      }
    }

    /* test for convergence */
    fsum = 0.0;
    for (j=0;j<=n;j++) {
      fsum += f[j];
    }
    favg = fsum/(n+1);
    s = 0.0;
    for (j=0;j<=n;j++) {
      s += pow((f[j]-favg),2.0)/(n);
    }
    s = sqrt(s);
    if (s < EPSILON) break;
  }
  /* end main loop of the minimization */

  /* find the index of the smallest value */
  vs=0;
  for (j=0;j<=n;j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  for (j = 0; j < n; j++){
    start[j] = v[vs][j];
  }
  //min = call_r_function(r_func, n1, p, v[vs], index, X, z_g, mu_g, sigma);

  // Clean up memory
  for(i = 0; i <= n; i++) delete[] v[i];
  delete[] v;
  delete[] f;
  delete[] vr;
  delete[] ve;
  delete[] vc;
  delete[] vm;
  free_carray(X, X_mat.nrow());
  free_carray(sigma, sigma_mat.nrow());

  return NumericVector(start_vec.begin(), start_vec.end());
}

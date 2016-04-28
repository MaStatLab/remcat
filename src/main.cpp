#include "remcat.h"
#include "optree.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

extern "C" {


  SEXP remcat_C(SEXP xg_mat, SEXP xe_mat, SEXP sample, SEXP kg, SEXP ke, SEXP rho0, SEXP mode);

};

SEXP remcat_C(SEXP xg_mat, SEXP xe_mat, SEXP sample, SEXP kg, SEXP ke, SEXP rho0, SEXP mode) {

  int n_PROTECTED = 0;

  PROTECT(xe_mat = coerceVector(xe_mat, INTSXP)); n_PROTECTED++;
  PROTECT(xg_mat = coerceVector(xg_mat, INTSXP)); n_PROTECTED++;
  PROTECT(sample = coerceVector(sample, INTSXP)); n_PROTECTED++;
  PROTECT(kg = coerceVector(kg, INTSXP)); n_PROTECTED++;
  PROTECT(ke = coerceVector(ke, INTSXP)); n_PROTECTED++;
  PROTECT(rho0 = AS_NUMERIC(rho0)); n_PROTECTED++;
  PROTECT(mode = coerceVector(mode,INTSXP)); n_PROTECTED++;

  COND_GBT my_cond_gbt(INTEGER(kg)[0],INTEGER(ke)[0],REAL(rho0)[0],INTEGER(mode)[0]);

  int n;
  int p;
  int i,j;

  int xg_curr_original;

  COV_TYPE xe_curr = 0;
  GENE_TYPE xg_curr;

  SEXP logrho;
  SEXP loggamma;
  SEXP logpsi;
  SEXP logphi;
  SEXP ans;


  PROTECT(ans = allocVector(VECSXP,4)); n_PROTECTED++;

  p = INTEGER(ke)[0];
  n = LENGTH(xe_mat)/p;

  for ( i = 0; i < n; i++ ) {
     xe_curr = 0;

    for ( j = 0; j < p; j++ ) {
      xe_curr = (xe_curr << 1) | (INTEGER(xe_mat)[n*j+i] & 1); // pth variable of ith individual
    }

    xg_curr.first = 0;
    xg_curr.second = 0;

    for ( j = 0; j < INTEGER(kg)[0]; j++ ) {
      xg_curr_original = INTEGER(xg_mat)[n*j+i];
      xg_curr.first = (xg_curr.first << 1) | ( xg_curr_original == 0); // pth variable of ith individual
      xg_curr.second = (xg_curr.second << 1) | (xg_curr_original == 2); // pth variable of ith individual

    }

    my_cond_gbt.add_data(xg_curr,xe_curr,INTEGER(sample)[i]);
  }


  my_cond_gbt.update();

  PROTECT(logphi = allocVector(REALSXP,1));n_PROTECTED++;
  PROTECT(logpsi = allocVector(REALSXP,1));n_PROTECTED++;
  PROTECT(logrho = allocVector(REALSXP,1));n_PROTECTED++;
  PROTECT(loggamma = allocVector(REALSXP,1));n_PROTECTED++;

  REAL(logphi)[0] = my_cond_gbt.get_root_logphi();
  REAL(logpsi)[0] = my_cond_gbt.get_root_logpsi();
  REAL(logrho)[0] = my_cond_gbt.get_root_logrho();
  REAL(loggamma)[0] = my_cond_gbt.get_root_loggamma();

  SET_VECTOR_ELT(ans, 0, logrho);
  SET_VECTOR_ELT(ans, 1, loggamma);
  SET_VECTOR_ELT(ans, 2, logphi);
  SET_VECTOR_ELT(ans, 3, logpsi);


  UNPROTECT(n_PROTECTED);


  return ans;

}



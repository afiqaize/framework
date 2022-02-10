// -*- C++ -*-
// author: afiq anuar
// short: codes relating to the usage of gsl chi2 fitter

#ifndef FWK_FIT_UTIL_H
#define FWK_FIT_UTIL_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_multifit.h"

class Workspace {
  gsl_multifit_linear_workspace *wsp;

public:
  Workspace(const size_t npoint, const size_t nparam): wsp(gsl_multifit_linear_alloc(npoint, nparam)) {}
  ~Workspace() { gsl_multifit_linear_free(wsp); }

  operator gsl_multifit_linear_workspace* () { return wsp; }
};



class Matrix {
  gsl_matrix *mat;

public:
  Matrix(const size_t nrow, const size_t ncol): mat(gsl_matrix_calloc(nrow, ncol)) {}
  ~Matrix() { gsl_matrix_free(mat); }

  double get(const size_t irow, const size_t icol) { return gsl_matrix_get(mat, irow, icol); }
  void set(const size_t irow, const size_t icol, double value) { gsl_matrix_set(mat, irow, icol, value); }

  operator gsl_matrix* () { return mat; }
};



class Vector {
  gsl_vector *vec;

public:
  Vector(const size_t ndim): vec(gsl_vector_calloc(ndim)) {}
  ~Vector() { gsl_vector_free(vec); }

  double get(const size_t idim) { return gsl_vector_get(vec, idim); }
  void set(const size_t idim, double value) { gsl_vector_set(vec, idim, value); }

  operator gsl_vector* () { return vec; }
};

#endif

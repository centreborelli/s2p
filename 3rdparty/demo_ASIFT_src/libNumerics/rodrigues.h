// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#ifndef RODRIGUES_H
#define RODRIGUES_H

#include "matrix.h"
#include <math.h>

namespace libNumerics {

/// Skew-symmetric matrix of 3-vector v.
template <class T> matrix<T> skew(const vector<T>& v);
/// Rodrigues's rotation: exp(w_x).
template <class T> matrix<T> rotation(vector<T> w);
/// Inverse Rodrigues's formula: w s.t. R=exp(w_x).
template <class T> vector<T> rotationAxis(const matrix<T>& R);

} // libNumerics

#include "rodrigues.cpp"

#endif

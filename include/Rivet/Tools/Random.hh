// -*- C++ -*-
#ifndef RIVET_Random_HH
#define RIVET_Random_HH

#include <random>
// #if defined(_OPENMP)
// #include "omp.h"
// #endif

namespace Rivet {

  /// Return a thread-safe random number generator (mainly for internal use)
  std::mt19937& rng();

  /// Return a uniformly sampled random number between 0 and 1
  double rand01();

  /// Return a random number sampled from a Gaussian/normal distribution
  double randnorm(double loc, double scale);

  /// Return a random number sampled from a log-normal distribution
  double randlognorm(double loc, double scale);

  /// Return a random number sampled from a Crystal Ball distribution
  double randcrystalball(double alpha, double n, double mu, double sigma);


  /// Probability density of a Gaussian/normal distribution at x
  double pNorm(double x, double mu, double sigma);
  /// Probability density of a Crystal Ball distribution at x
  double pCrystalBall(double x, double alpha, double n, double mu, double sigma);


}

#endif

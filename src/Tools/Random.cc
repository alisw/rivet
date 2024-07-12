// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include <random>
#if defined(_OPENMP)
#include "omp.h"
#endif

namespace Rivet {

  using namespace std;


  // Return a thread-safe random number generator
  mt19937& rng() {
    #if defined(_OPENMP)
    static map<int,mt19937> gens;
    const int nthread = omp_get_thread_num();
    if (gens.find(nthread) == gens.end()) {
      // Make seeds for each thread, either via the standard seed generator or based on a fixed seed from the environment
      vector<uint32_t> seeds(nthread+1);
      const uint32_t envseed = getEnvParam<uint32_t>("RIVET_RANDOM_SEED", 0);
      //cout << "RIVET_RANDOM_SEED = " << envseed << endl;
      if (envseed > 0) {
        std::iota(seeds.begin(), seeds.end(), envseed);
      } else {
        seed_seq seq{1,2,3,4,5};
        seq.generate(seeds.begin(), seeds.end());
      }
      gens[nthread] = mt19937(seeds[nthread]);
      //cout << "Thread " << nthread+1 << ", seed=" << seeds[nthread] << " (" << gens.size() << " RNGs)" << endl;
    }
    mt19937& g = gens[nthread];
    #else
    const uint32_t envseed = getEnvParam<uint32_t>("RIVET_RANDOM_SEED", 12345);
    static mt19937 g(envseed);
    #endif
    return g;
  }


  // Return a uniformly sampled random number between 0 and 1
  double rand01() {
    const double x = generate_canonical<double, 32>(rng()); ///< @todo What's the "correct" number of bits of randomness?
    //cout << "RAND01 -> " << x << endl;
    return x;
  }


  // Return a Gaussian/normal sampled random number with the given mean and width
  double randnorm(double loc, double scale) {
    normal_distribution<> d(loc, scale);
    const double x = d(rng());
    //cout << "RANDNORM -> " << x << endl;
    return x;
  }


  // Return a log-normal sampled random number
  double randlognorm(double loc, double scale) {
    lognormal_distribution<> d(loc, scale);
    const double x = d(rng());
    //cout << "RANDLOGNORM -> " << x << endl;
    return x;
  }


  // Return a random number sampled from a Crystal Ball distribution
  double randcrystalball(double alpha, double n, double mu, double sigma) {
    const double aalpha = fabs(alpha);
    const double nalpha = n/aalpha;
    const double C = nalpha / (n-1) * exp(-sqr(aalpha)/2.);
    const double D = sqrt(M_PI/2) * (1 + std::erf(aalpha/M_SQRT2));
    const double normfrac = D / (C+D);
    if (rand01() < normfrac) { // Sample from (the relevant bit of) the Gaussian
      while (1) { // sample until we get a value that isn't in the power-law tail
        const double x = randnorm(mu, sigma);
        if (x-mu >= -alpha*sigma) return x;
      }
    } else { // Sample from the power-law tail
      // Cf. https://stackoverflow.com/questions/918736/random-number-generator-that-produces-a-power-law-distribution
      const double xt = nalpha * pow(1-rand01(), 1/(1-n));
      const double x = sigma*(nalpha - xt - aalpha) + mu;
      return x;
    }
    return NAN;
  }



  double pNorm(double x, double mu, double sigma) {
    const double dx = x - mu;
    const double y = dx/sigma;
    const double p = exp(-y*y/2.);
    return p / sqrt(TWOPI) / sigma;
  }

  double pCrystalBall(double x, double alpha, double n, double mu, double sigma) {
    const double dx = x - mu;
    const double y = dx/sigma;
    const double aalpha = fabs(alpha);
    const double nalpha = n/aalpha;
    double p = -1;
    if (y < -alpha) { // CB part
      const double A = pow(nalpha,n) * exp(-sqr(aalpha)/2.);
      const double B = nalpha - aalpha;
      p = A * pow(B-y, -n);
    } else { // Gaussian part
      p = exp(-y*y/2.);
    }
    // Normalize
    const double C = nalpha / (n-1) * exp(-sqr(aalpha)/2.);
    const double D = sqrt(M_PI/2) * (1 + std::erf(aalpha/M_SQRT2));
    const double Z = sigma * (C + D);
    return p/Z;
  }


}

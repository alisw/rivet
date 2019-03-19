// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include <random>
#if defined(_OPENMP)
#include "omp.h"
#endif

namespace Rivet {


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
    static mt19937 g(12345);
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


}

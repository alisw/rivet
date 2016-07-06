// -*- C++ -*-
#ifndef RIVET_MT2_HH
#define RIVET_MT2_HH

namespace Rivet {

  class FourMomentum;

  namespace mT2 {

    /// @brief Convenience wrapper for the mT2 calculator of Cheng/Han
    ///
    /// (arXiv:0810.5178). Could be adapted for other backends in future.
    double mT2(const FourMomentum& a,
               const FourMomentum& b,
               const FourMomentum& pTmiss,
               double invisiblesMass);

  }

}

#endif

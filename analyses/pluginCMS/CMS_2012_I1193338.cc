// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2012_I1193338 : public Analysis {
  public:

    CMS_2012_I1193338()
      : Analysis("CMS_2012_I1193338")
    {    }


    void init() {
      declare(ChargedFinalState((Cuts::etaIn(-2.4, 2.4) && Cuts::pT >=  0.2*GeV)), "CFS");
      declare(FinalState(), "FS");

      book(_h_sigma ,1, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = 1.0;

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() > 1) {_h_sigma->fill(1.5, weight);}
      if (cfs.size() > 2) {_h_sigma->fill(2.5, weight);}
      if (cfs.size() > 3) {_h_sigma->fill(3.5, weight);}

      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.size() < 2) vetoEvent; // need at least two particles to calculate gaps

      double gapcenter = 0.;
      double LRG = 0.;
      double etapre = 0.;
      bool first = true;

      for(const Particle& p : fs.particles(cmpMomByEta)) { // sorted from minus to plus
        if (first) { // First particle
          first = false;
          etapre = p.eta();
        } else {
          double gap = fabs(p.eta()-etapre);
          if (gap > LRG) {
            LRG = gap; // largest gap
            gapcenter = (p.eta()+etapre)/2.; // find the center of the gap to separate the X and Y systems.
          }
          etapre = p.eta();
        }
      }


      FourMomentum mxFourVector, myFourVector;
      for(const Particle& p : fs.particles(cmpMomByEta)) {
        ((p.eta() > gapcenter) ? mxFourVector : myFourVector) += p.momentum();
      }
      const double M2 = max(mxFourVector.mass2(), myFourVector.mass2());
      const double xi = M2/sqr(sqrtS()); // sqrt(s)=7000 GeV, note that units cancel
      if (xi < 5e-6) vetoEvent;

      _h_sigma->fill(0.5, weight);
    }


    void finalize() {
      scale(_h_sigma, crossSection()/millibarn/sumOfWeights());
    }

  private:

    Histo1DPtr _h_sigma;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1193338);

}

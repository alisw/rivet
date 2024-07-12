// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  /// Production cross-sections of muons from $b$ hadron decays in $pp$ collisions
  class CMS_2011_S8941262 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S8941262);


    /// Book histograms and initialise projections before the run
    void init() {
      book(_h_total ,1, 1, 1);
      book(_h_mupt  ,2, 1, 1);
      book(_h_mueta ,3, 1, 1);
      nbtot=0.;   nbmutot=0.;

      IdentifiedFinalState ifs(Cuts::abseta < 2.1 && Cuts::pT > 6*GeV);
      ifs.acceptIdPair(PID::MUON);
      declare(ifs, "IFS");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // a b-quark must have been produced
      /// @todo Ouch. Use hadron tagging...
      int nb = 0;
      for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())) {
        if (abs(p->pdg_id()) == PID::BQUARK) nb += 1;
      }
      if (nb == 0) vetoEvent;
      nbtot += weight;

      // Event must contain a muon
      Particles muons = apply<IdentifiedFinalState>(event, "IFS").particlesByPt();
      if (muons.size() < 1) vetoEvent;
      nbmutot += weight;

      FourMomentum pmu = muons[0].momentum();
      _h_total->fill(      7000/GeV, weight);
      _h_mupt->fill(   pmu.pT()/GeV, weight);
      _h_mueta->fill( pmu.eta()/GeV, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_total, crossSection()/microbarn/sumOfWeights());
      scale(_h_mupt,  crossSection()/nanobarn/sumOfWeights());
      scale(_h_mueta, crossSection()/nanobarn/sumOfWeights());
    }


  private:

    /// @todo Convert to counters?
    double nbtot, nbmutot;

    /// @{
    Histo1DPtr _h_total;
    Histo1DPtr _h_mupt;
    Histo1DPtr _h_mueta;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S8941262, CMS_2011_I884811);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"

namespace Rivet {


  /// @brief Just measures a few observables as a demo
  class EXAMPLE : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE);


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Projections
      const FinalState cnfs(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV);
      const ChargedFinalState cfs(cnfs);
      declare(cnfs, "FS");
      declare(cfs, "CFS");
      declare(FastJets(cnfs, FastJets::ANTIKT, 0.4), "Jets");
      declare(Thrust(cfs), "Thrust");
      declare(Sphericity(cfs), "Sphericity");

      // Histograms
      book(_histTot         ,"TotalMult", 100, -0.5, 99.5);
      book(_histChTot       ,"TotalChMult", 50, -1.0, 99.0);
      book(_histHadrTot     ,"HadrTotalMult", 100, -0.5, 99.5);
      book(_histHadrChTot   ,"HadrTotalChMult", 50, -1.0, 99.0);
      book(_histMajor       ,"Major", 10, 0.0, 0.6);
      book(_histSphericity  ,"Sphericity", 10, 0.0, 0.8);
      book(_histAplanarity  ,"Aplanarity", 10, 0.0, 0.3);
      book(_histThrust      ,"Thrust", { 0.5, 0.6, 0.7, 0.80, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0 });
    }


    /// Do the analysis
    void analyze(const Event& event) {
      const Particles& cnparticles = apply<FinalState>(event, "FS").particles();
      MSG_DEBUG("Total multiplicity = " << cnparticles.size());
      _histTot->fill(cnparticles.size());
      int cnhadronmult = 0;
      for (const Particle& p : cnparticles)
        if (isHadron(p)) cnhadronmult += 1;
      MSG_DEBUG("Hadron multiplicity = " << cnhadronmult);
      _histHadrTot->fill(cnhadronmult);

      const Particles& cparticles = apply<FinalState>(event, "CFS").particles();
      MSG_DEBUG("Total charged multiplicity = " << cparticles.size());
      _histChTot->fill(cparticles.size());
      int chadronmult = 0;
      for (const Particle& p : cparticles)
        if (isHadron(p)) chadronmult += 1;
      MSG_DEBUG("Hadron charged multiplicity = " << chadronmult);
      _histHadrChTot->fill(chadronmult);

      const Thrust& t = apply<Thrust>(event, "Thrust");
      MSG_DEBUG("Thrust = " << t.thrust());
      _histThrust->fill(t.thrust());
      _histMajor->fill(t.thrustMajor());

      const Sphericity& s = apply<Sphericity>(event, "Sphericity");
      MSG_DEBUG("Sphericity = " << s.sphericity());
      _histSphericity->fill(s.sphericity());
      MSG_DEBUG("Aplanarity = " << s.aplanarity());
      _histAplanarity->fill(s.aplanarity());

      const Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::pT > 5*GeV);
      const size_t num_b_jets = count_if(jets.begin(), jets.end(), hasBTag(Cuts::pT > 500*MeV));
      MSG_DEBUG("Num B-jets with pT > 5 GeV = " << num_b_jets);
    }


    /// Finalize
    void finalize() {
      normalize(_histTot); normalize(_histChTot); normalize(_histHadrTot); 
      normalize(_histHadrChTot); normalize(_histThrust); normalize(_histMajor); 
      normalize(_histSphericity); normalize(_histAplanarity);
    }

    //@}


    //@{
    /// Histograms
    Histo1DPtr _histTot, _histChTot, _histHadrTot, _histHadrChTot, _histThrust, _histMajor, _histSphericity, _histAplanarity;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(EXAMPLE);

}

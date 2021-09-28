// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  // This analysis is a derived from the class Analysis:
  class CMS_2013_I1208923 : public Analysis {

  public:
    // Constructor
    CMS_2013_I1208923()
      : Analysis("CMS_2013_I1208923") {

    }

    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;
      declare(fs, "FS");

      // Initialize the projections
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "Jets");

      // Book histograms
      {Histo1DPtr tmp; _h_sigma.add(0.0, 0.5,   book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_sigma.add(0.5, 1.0,   book(tmp, 1, 1, 2));}
      {Histo1DPtr tmp; _h_sigma.add(1.0, 1.5,   book(tmp, 1, 1, 3));}
      {Histo1DPtr tmp; _h_sigma.add(1.5, 2.0,   book(tmp, 1, 1, 4));}
      {Histo1DPtr tmp; _h_sigma.add(2.0, 2.5,   book(tmp, 1, 1, 5));}
      
      {Histo1DPtr tmp; _h_invMass.add(0.0, 0.5, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_invMass.add(0.5, 1.0, book(tmp, 2, 1, 2));}
      {Histo1DPtr tmp; _h_invMass.add(1.0, 1.5, book(tmp, 2, 1, 3));}
      {Histo1DPtr tmp; _h_invMass.add(1.5, 2.0, book(tmp, 2, 1, 4));}
      {Histo1DPtr tmp; _h_invMass.add(2.0, 2.5, book(tmp, 2, 1, 5));}
    }

    // Analysis
    void analyze(const Event &event) {
      const double weight = 1.0;
      const FastJets &fJets = apply<FastJets>(event, "Jets");
      
      // Fill the jet pT spectra
      const Jets& jets = fJets.jetsByPt(Cuts::pt>100.*GeV && Cuts::absrap <2.5);
      for (const Jet &j : jets) {
        _h_sigma.fill(fabs(j.momentum().rapidity()), j.momentum().pT() / GeV, weight);
      }

      // Require two jets
      const Jets& dijets = fJets.jetsByPt(Cuts::pt>30.*GeV && Cuts::absrap < 2.5);
      if (dijets.size() > 1) {
        if (dijets[0].momentum().pT() / GeV > 60.) {
          // Fill the invariant mass histogram
          double ymax = max(dijets[0].momentum().absrapidity(), dijets[1].momentum().absrapidity());
          double invMass = FourMomentum(dijets[0].momentum() + dijets[1].momentum()).mass();
          _h_invMass.fill(fabs(ymax), invMass, weight);
        }
      } 

    }


    // Scale histograms by the production cross section
    void finalize() {
      _h_sigma.scale(  crossSection() / sumOfWeights() / 2.0, this);
      _h_invMass.scale(crossSection() / sumOfWeights() / 2.0, this);
    }

  private:
    BinnedHistogram _h_sigma;
    BinnedHistogram _h_invMass;
  };

  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2013_I1208923);
}

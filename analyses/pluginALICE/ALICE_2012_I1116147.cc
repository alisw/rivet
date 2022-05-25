//-*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class ALICE_2012_I1116147 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2012_I1116147);


    /// Initialise projections and histograms
    void init() {

      const UnstableParticles ufs(Cuts::absrap < RAPMAX);
      declare(ufs, "UFS");

      // Check if cm energy is 7 TeV or 0.9 TeV
      if (isCompatibleWithSqrtS(900))       _cm_energy_case = 1;
      else if (isCompatibleWithSqrtS(7000)) _cm_energy_case = 2;
      if (_cm_energy_case == 0)
        throw UserError("Center of mass energy of the given input is neither 900 nor 7000 GeV.");

      // Book histos
      if (_cm_energy_case == 1) {
        book(_h_pi0,       2,1,1);
      } else {
        book(_h_pi0,       1,1,1);
        book(_h_eta,       3,1,1);
        book(_h_etaToPion, 4,1,1);
      }

      // Temporary plots with the binning of _h_etaToPion to construct the eta/pi0 ratio
      book(_temp_h_pion, "TMP/h_pion", refData(4,1,1));
      book(_temp_h_eta , "TMP/h_eta",  refData(4,1,1));
    }


    /// Per-event analysis
    void analyze(const Event& event) {

      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const double normfactor = TWOPI*p.pT()/GeV*2*RAPMAX;
        if (p.pid() == 111) {
          // Neutral pion; ALICE corrects for pi0 feed-down from K_0_s and Lambda
          if (p.hasAncestor(310) || p.hasAncestor(3122) || p.hasAncestor(-3122)) continue; //< K_0_s, Lambda, Anti-Lambda
          _h_pi0->fill(p.pT()/GeV, 1.0/normfactor);
          _temp_h_pion->fill(p.pT()/GeV);
        } else if (p.pid() == 221 && _cm_energy_case == 2) {
          // eta meson (only for 7 TeV)
          _h_eta->fill(p.pT()/GeV, 1.0/normfactor);
          _temp_h_eta->fill(p.pT()/GeV);
        }
      }
    }


    /// Normalize histos and construct ratio
    void finalize() {
      scale(_h_pi0, crossSection()/microbarn/sumOfWeights());
      if (_cm_energy_case == 2) {
        divide(_temp_h_eta, _temp_h_pion, _h_etaToPion);
        scale(_h_eta, crossSection()/microbarn/sumOfWeights());
      }
    }


  private:

    const double RAPMAX = 0.8;
    int _cm_energy_case = 0;

    Histo1DPtr _h_pi0, _h_eta;
    Histo1DPtr _temp_h_pion, _temp_h_eta;
    Scatter2DPtr _h_etaToPion;

  };


  RIVET_DECLARE_PLUGIN(ALICE_2012_I1116147);

}

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Z pT and Z phi* at 13 TeV
  class ATLAS_2019_I1768911 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2019_I1768911);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options 
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      if ( getOption("LMODE") == "MU" ) _mode = 2;

      // Configure projections
      const FinalState fs;
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 27*GeV;
      ZFinder zmmFinder(fs, cuts, PID::MUON, 66*GeV, 116*GeV, 0.1,
		      ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
      declare(zmmFinder, "ZFinder_mu");
      ZFinder zeeFinder(fs, cuts, PID::ELECTRON, 66*GeV, 116*GeV, 0.1,
		      ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
      declare(zeeFinder, "ZFinder_el");

      // Book histograms
      book(_h["zpt_combined_dressed_normalised"], 27, 1, 1);
      book(_h["zphistar_combined_dressed_normalised"], 28, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get leptonic Z boson
      const ZFinder& zmmFinder = apply<ZFinder>(event, "ZFinder_mu");
      const ZFinder& zeeFinder = apply<ZFinder>(event, "ZFinder_el");
      if (_mode == 2 && zmmFinder.bosons().size() != 1 && zeeFinder.bosons().size())   vetoEvent;
      if (_mode == 1 && zeeFinder.bosons().size() != 1 && zmmFinder.bosons().size())   vetoEvent;
      if (_mode == 0 && (zeeFinder.bosons().size() + zmmFinder.bosons().size()) != 1)  vetoEvent;
      const Particle& Zboson = zeeFinder.bosons().size()? zeeFinder.boson() : zmmFinder.boson();

      // cut on Z boson leptons and calculate pTll and phistar
      const Particles& leptons = zeeFinder.bosons().size()? zeeFinder.constituents() : zmmFinder.constituents();
      if (leptons.size() != 2 || leptons[0].charge3() * leptons[1].charge3() > 0) vetoEvent;

      const double zpt = Zboson.pT()/GeV;
      const Particle& lminus = leptons[0].charge() < 0 ? leptons[0] : leptons[1];
      const Particle& lplus  = leptons[0].charge() < 0 ? leptons[1] : leptons[0];
      const double phi_acop = M_PI - deltaPhi(lminus, lplus);
      const double costhetastar = tanh( 0.5 * (lminus.eta() - lplus.eta()) );
      const double sin2thetastar = (costhetastar > 1) ? 0.0 : (1.0 - sqr(costhetastar));
      const double phistar = tan(0.5 * phi_acop) * sqrt(sin2thetastar);

      _h["zpt_combined_dressed_normalised"]->fill(zpt);
      _h["zphistar_combined_dressed_normalised"]->fill(phistar);
    }


    /// Normalise histograms etc., after the run
    void finalize() {  normalize(_h);  }
    //@}


  protected:

    size_t _mode;


  private:

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    //@}

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2019_I1768911);
}

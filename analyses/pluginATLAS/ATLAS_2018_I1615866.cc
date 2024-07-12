// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> mu+ mu-
  class ATLAS_2018_I1615866 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1615866);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      ChargedFinalState cfs(Cuts::abseta < 2.5 and Cuts::pT>0.4);
      declare(cfs,"CFS");
      // Get muons which pass the initial kinematic cuts
      IdentifiedFinalState muon_fs(Cuts::abseta < 2.4 && Cuts::pT > 6*GeV);
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "MUON_FS");
      // histograms
      book(_h_sigma,1,1,1);
      book(_h_mass ,2,1,1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() != 2) vetoEvent; // no other charged particles in 2.5
      // extract the muons
      const Particles& muonFS = apply<IdentifiedFinalState>(event, "MUON_FS").particles();
      // have 2 muons with opposite charge
      if (muonFS.size() != 2) vetoEvent;
      if (muonFS[0].pid() != -muonFS[1].pid()) vetoEvent;
      // invariant mass between 12 and 70
      double mmumu = (muonFS[0].momentum()+muonFS[1].momentum()).mass();
      if(mmumu<12. or mmumu>70.) vetoEvent;
      // cut pt >10 if pair mass >30
      for(unsigned int ix=0;ix<2;++ix) {
	if(mmumu>30 && muonFS[ix].perp()<10.) vetoEvent;
      }
      // fill histogram
      _h_mass->fill(mmumu);
      _h_sigma->fill(13000);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sigma, crossSection()/picobarn/sumOfWeights());
      scale(_h_mass , crossSection()/picobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma;
    Histo1DPtr _h_mass;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1615866);

}

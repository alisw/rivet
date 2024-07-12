// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> l+ l-
  class ATLAS_2015_I1377585 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2015_I1377585);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      ChargedFinalState cfs(Cuts::abseta < 2.5 and Cuts::pT>0.4);
      declare(cfs,"CFS");
      // Get electrons which pass the initial kinematic cuts
      IdentifiedFinalState electron_fs(Cuts::abseta < 2.4 && Cuts::pT > 12*GeV);
      electron_fs.acceptIdPair(PID::ELECTRON);
      declare(electron_fs, "ELECTRON_FS");
      // Get muons which pass the initial kinematic cuts
      IdentifiedFinalState muon_fs(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "MUON_FS");
      // book histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_sigma [ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() != 2) vetoEvent; // no other charged particles in 2.5
      // have e+e-
      const Particles& electronFS = apply<IdentifiedFinalState>(event, "ELECTRON_FS").particles();
      if (electronFS.size() == 2 && electronFS[0].pid() == -electronFS[1].pid()) {
	double mee = (electronFS[0].momentum()+electronFS[1].momentum()).mass();
	if(mee>24.) _h_sigma[0]->fill(7000);
      }
      // have 2 muons with opposite charge
      const Particles& muonFS = apply<IdentifiedFinalState>(event, "MUON_FS").particles();
      if (muonFS.size() == 2 && muonFS[0].pid() == -muonFS[1].pid()) {
	double mmumu = (muonFS[0].momentum()+muonFS[1].momentum()).mass();
	if(mmumu>20.) _h_sigma[1]->fill(7000);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	scale(_h_sigma [ix], crossSection()/picobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2015_I1377585);

}

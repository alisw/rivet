// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// Exclusive photon-photon production of muon pairs at 7 TeV
  class CMS_2011_I954992 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2011_I954992);


    void init() {
      ChargedFinalState cfs(Cuts::abseta < 2.4);
      declare(cfs,"CFS");

      // Get muons which pass the initial kinematic cuts
      IdentifiedFinalState muon_fs(Cuts::abseta < 2.1 && Cuts::pT > 4*GeV);
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "MUON_FS");

      _h_sigma = bookHisto1D(1,1,1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() != 2) vetoEvent; // no other charged particles in 2.4

      const Particles& muonFS = apply<IdentifiedFinalState>(event, "MUON_FS").particles();
      if (muonFS.size() != 2) vetoEvent;

      if (muonFS[0].charge() != muonFS[1].charge()) {
         const double dimuon_mass = (muonFS[0].momentum() + muonFS[1].momentum()).mass();
         const double v_angle     = muonFS[0].momentum().angle(muonFS[1].momentum());
         const double dPhi        = deltaPhi(muonFS[0], muonFS[1]);
         const double deltaPt     = fabs(muonFS[0].pT() - muonFS[1].pT());

         if (dimuon_mass >= 11.5*GeV &&
             v_angle < 0.95*PI       &&
             dPhi    > 0.9*PI        &&
             deltaPt < 1.*GeV        ) {
           _h_sigma->fill(sqrtS()/GeV, weight);
         }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sigma, crossSection()/picobarn/sumOfWeights());
    }


    /// Histogram
    Histo1DPtr _h_sigma;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_I954992);

}

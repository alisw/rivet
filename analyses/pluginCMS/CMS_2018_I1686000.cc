// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// Fiducial single-top + photon cross-section measurement at 13 TeV
  class CMS_2018_I1686000 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1686000);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Leptons
      declare(DressedLeptons(PromptFinalState(), 0.1, Cuts::abseta < 2.4 && Cuts::pT > 26*GeV), "Leptons");

      // Jets
      declare(FastJets(FinalState(Cuts::abseta < 5), FastJets::ANTIKT, 0.4), "Jets");

      // Photons
      declare(PromptFinalState(Cuts::pid == PID::PHOTON && Cuts::pT > 25*GeV && Cuts::abseta < 1.44), "Photons");

      // MET
      declare(MissingMomentum(FinalState(Cuts::abseta < 5)), "MET");


      // Book xsec counter
      book(_c_xsec_fid, "xsec");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // // Find at least 2 jets, one b-tagged
      // const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::abseta < 4.7 && Cuts::pT > 40*GeV);
      // Jets bjets, ljets;
      // for (const Jet& j : jets)
      //   ((j.abseta() < 2.5 && j.bTagged()) ? bjets : ljets) += j;
      // if (bjets.empty() || ljets.empty()) vetoEvent;
      // const Jet& bjet = bjets[0];
      // const Jet& ljet = ljets[0];

      // // Require exactly one isolated lepton, and it has to be a muon
      // const Particles& leps = apply<FinalState>(event, "Leptons").particlesByPt();
      // const Particles isoleps = discardIfAnyDeltaRLess(leps, jets, 0.3);
      // if (isoleps.size() != 1 || isoleps[0].abspid() != PID::MUON) vetoEvent;
      // const Particle& muon = isoleps[0];

      // // Require exactly one isolated photon
      // const Particles& photons = apply<FinalState>(event, "Photons").particlesByPt();
      // const Particles muisophotons = filter_discard(photons, deltaRLess(muon,0.5));
      // const Particles isophotons = discardIfAnyDeltaRLess(muisophotons, Jets{bjet,ljet}, 0.5);
      // if (isophotons.size() != 1) vetoEvent;

      // // Require 30 GeV of missing ET
      // const double met = apply<MissingMomentum>(event, "MET").met();
      // if (met < 30*GeV) vetoEvent;


      // Find light jets
      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt();
      const Jets ljets = filter_discard(jets, [](const Jet& j){ return j.abseta() < 2.5 && j.bTagged(); } );

      // Require a photon, isolated from the light jet
      Particles photons = apply<FinalState>(event, "Photons").particlesByPt();
      if (!ljets.empty()) ifilter_discard(photons, deltaRLess(ljets[0], 0.5));
      if (photons.empty()) vetoEvent;

      // Fill counter
      _c_xsec_fid->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double BRmu = 0.13 + 0.13*0.17; //< decay BR direct to a muon, and to a muon via a tau
      scale(_c_xsec_fid, BRmu*crossSection()/femtobarn/sumOfWeights());
    }

    //@}


    /// Counter
    CounterPtr _c_xsec_fid;


  };


  RIVET_DECLARE_PLUGIN(CMS_2018_I1686000);

}

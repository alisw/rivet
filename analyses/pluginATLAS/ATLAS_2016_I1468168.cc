#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  class ATLAS_2016_I1468168 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1468168);

    void init() {
      // Eta ranges
      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT >= 1.0*MeV;
      // Lepton cuts
      Cut lep_cuts = Cuts::abseta < 2.5 && Cuts::pT >= 25.0*GeV;

      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts, true);
      declare(dressedelectrons, "DressedElectrons");

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true);
      declare(dressedmuons, "DressedMuons");

      book(_s , 2, 1 ,1);
      book(_c , "_counter");
    }


    void analyze(const Event& event) {

      // Get the selected objects, using the projections.
      const size_t num_es = apply<DressedLeptons>(event, "DressedElectrons").dressedLeptons().size();
      const size_t num_mus = apply<DressedLeptons>(event, "DressedMuons").dressedLeptons().size();

      // Evaluate basic event selection
      const bool pass_emu = num_es == 1 && num_mus == 1;
      if (!pass_emu) vetoEvent;

      // Fill histogram to measure the event acceptance
      _c->fill();
    }


    void finalize() {
      // Normalize to cross-section
      scale(_c, crossSection() / sumOfWeights());

      double err = _c->err();
      _s->addPoint(13*TeV, _c->val(), make_pair(0.5, 0.5), make_pair(err,err));

    }


   private:

    CounterPtr _c;
    Scatter2DPtr _s;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1468168);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {


  /// @brief Colour reconnection in ttbar dilepton events
  class ATLAS_2022_I2152933 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2022_I2152933);

    void book2D(std::string name, std::vector<double> doubleDiff_bins, unsigned int table){
    	for (unsigned int i = 0; i < doubleDiff_bins.size() - 1; ++i){
				std::string nbin = std::to_string(i);
				std::string title = name+"_"+nbin;
				Histo1DPtr tmp;
				_h_multi[name].add(doubleDiff_bins[i], doubleDiff_bins[i+1], book(tmp, table+i,1,1));
			}
    }

    double integral2D(BinnedHistogram& h_multi) {
        double total_integral = 0;
        for  (Histo1DPtr& h : h_multi.histos()){
          total_integral += h->integral(false);
        }
        return total_integral;
    }


    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

    	// Photons
    	PromptFinalState photons(Cuts::abspid == PID::PHOTON);

    	// Muons
    	Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
    	PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true);
    	DressedLeptons all_dressed_mu(photons, bare_mu, 0.1, lepton_cuts, true);
    	declare(all_dressed_mu, "muons");

    	// Electrons
    	PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true);
    	DressedLeptons all_dressed_el(photons, bare_el, 0.1, lepton_cuts, true);
    	declare(all_dressed_el, "electrons");

    	//Jet forming
    	const InvisibleFinalState neutrinos(true, true);

    	VetoedFinalState vfs(FinalState(Cuts::abseta < 5.0));
    	vfs.addVetoOnThisFinalState(all_dressed_el);
    	vfs.addVetoOnThisFinalState(all_dressed_mu);
    	vfs.addVetoOnThisFinalState(neutrinos);
    	FastJets jetfs(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
    	declare(jetfs, "jets");

      // FinalState charged particles
      ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      declare(tracks, "tracks");

      declare(MissingMomentum(), "ETmiss");

      // Book histograms
      book(_h["nch"], 7,1,1);
      book(_h["sumPt"], 8,1,1);

      book2D("sumPt_nch_multi", nch_2D_bins, 9);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "muons").dressedLeptons();
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      // Retrieve charge final state particles, sorted by pT
      const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();

      // OVERLAP REMOVAL
      idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      idiscardIfAnyDeltaRLess(electrons, jets, 0.4);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV);
      });

      // Veto event if there are not exactly one electron and one muon
      if ( electrons.size() != 1 || muons.size() != 1)  vetoEvent;

      // Veto event if the selected electron and muon are not OS
      if (electrons[0].charge() == muons[0].charge())  vetoEvent;

      // Veto event if there are no 2 jets or 3 jets
      if (jets.size() < 2 || jets.size() > 3)  vetoEvent;

      // Veto event if the number of b-jets not exactly 2
      if (bjets.size() != 2)  vetoEvent;

      // Veto event with dilepton mass < 15 GeV
      FourMomentum ll = electrons[0].momentum() + muons[0].momentum();
      if (ll.mass() <= 15*GeV ) vetoEvent;

      // Remove tracks inside jets and leptons: within 0.4 delatR from a jet and 0.01 from a lepton

      Particles myTracks;
      for (const Particle& p : tracks) {
        bool isJetAssoc = false;
        for (unsigned int i = 0; i < jets.size(); ++i) {
          if(deltaR(jets[i], p, PSEUDORAPIDITY) < 0.4) {
          	isJetAssoc= true;
          	break;
          }
        }
        if (isJetAssoc) continue;
        if(deltaR(muons[0], p, PSEUDORAPIDITY) < 0.01) continue;
        if(deltaR(electrons[0], p, PSEUDORAPIDITY) < 0.01) continue;
        if(p.pid() == 0 || p.isStable() != 1 || p.charge() == 0) continue;
        myTracks.push_back(p);
      }

      // Veto event if the number of selected charged particles higher than 100
      if (myTracks.size() > 100)  vetoEvent;

      // Fill nch, sumpt, and sumpt vs. nch
      const double nch = min(myTracks.size(), 99u); // put overflow in the last bin
      _h["nch"]->fill(nch);

      const double sumPt = min(sum(myTracks, Kin::pT, 0.0), 119*GeV);
      _h["sumPt"]->fill(sumPt/GeV);

      const double sumPt2d = min(sumPt, 99*GeV);
      _h_multi["sumPt_nch_multi"].fill(nch, sumPt2d/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Normalize to unity
      normalize(_h);
      for (auto& hist : _h_multi) {
      	// scaling for normalised distribution according integral of whole set
        const double norm2D = integral2D(hist.second);
        hist.second.scale(1./norm2D, this);
      }
    }


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    map<string, BinnedHistogram> _h_multi;

    const vector<double> nch_2D_bins = { 0, 19.5, 39.5, 59.5, 79.5, 101 };
    /// @}

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2022_I2152933);

}

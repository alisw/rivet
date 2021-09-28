// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

  /// @brief Z + b(b) in pp at 13 TeV
  class ATLAS_2020_I1788444 : public Analysis {
  public:
    /// @name Constructors etc.
    //@{
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2020_I1788444);
    //@}

    /// Book histograms and initialise projections before the run
    void init() {

      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      if ( getOption("LMODE") == "MU" ) _mode = 2;

	    const FinalState fs;
	    // Define fiducial cuts for the leptons in the ZFinder
	    Cut lepcuts = (Cuts::pT > 27*GeV) & (Cuts::abseta < 2.5);
	    ZFinder zfinderE(fs, lepcuts, PID::ELECTRON, 76*GeV, 106*GeV);
	    ZFinder zfinderM(fs, lepcuts, PID::MUON, 76*GeV, 106*GeV);
	    declare(zfinderE, "zfinderE");
      declare(zfinderM, "zfinderM");

	    declare(HeavyHadrons(), "HFHadrons");

	    // // Photons
	    FinalState photons(Cuts::abspid == PID::PHOTON);
	    // Muons
	    PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true);
	    DressedLeptons all_dressed_mu(photons, bare_mu, 0.1, Cuts::abseta < 2.5, true);
	    // Electrons
	    PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true);
	    DressedLeptons all_dressed_el(photons, bare_el, 0.1, Cuts::abseta < 2.5, true);

	    //Jet forming
	    VetoedFinalState vfs(FinalState(Cuts::abseta < 4.5));
	    vfs.addVetoOnThisFinalState(all_dressed_el);
	    vfs.addVetoOnThisFinalState(all_dressed_mu);

	    FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
	    declare(jets, "jets");

	    // book histos - binning taken from data.yoda
      book(_h["i1b_ZpT"],2,1,1);
      book(_h["i1b_ZY"],4,1,1);
      book(_h["i1b_dPhiZb"],6,1,1);
      book(_h["i1b_dRZb"],8,1,1);
      book(_h["i1b_dYZb"],7,1,1);
      book(_h["i1b_bpT"],3,1,1);
      book(_h["i1b_bY"],5,1,1);

      book(_h["i2b_ZpT"],13,1,1);
      book(_h["i2b_dPhibb"],9,1,1);
      book(_h["i2b_dRbb"],11,1,1);
      book(_h["i2b_dYbb"],10,1,1);
      book(_h["i2b_Mbb"],12,1,1);
      book(_h["i2b_pTbb"],14,1,1);
      book(_h["i2b_pTOnMbb"],15,1,1);

      book(_h["ib_nBJets"],1,1,1);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinderE = apply<ZFinder>(event, "zfinderE");
      const Particles& els = zfinderE.constituents();
	    const ZFinder& zfinderM = apply<ZFinder>(event, "zfinderM");
	    const Particles& mus = zfinderM.constituents();

      // default is to run average of Z->ee and Z->mm
      // use LMODE option to pick one channel
	    if ( (els.size() + mus.size()) != 2 )  vetoEvent;

      if (      _mode == 0 && !(els.size()==2 || mus.size()==2) )  vetoEvent;
	    else if ( _mode == 1 && !(els.size() == 2 && mus.empty()) )  vetoEvent;
	    else if ( _mode == 2 && !(els.empty() && mus.size() == 2) )  vetoEvent;

	    double Vpt = 0, Vy = 0, Veta = 0, Vphi = 0;

	    if ( els.size()==2 ) {
        Vpt = zfinderE.boson().pt()/GeV;
        Vphi = zfinderE.boson().phi();
        Vy = zfinderE.boson().rapidity();
        Veta = zfinderE.boson().eta();
	    } else {
        Vpt = zfinderM.boson().pt()/GeV;
        Vphi = zfinderM.boson().phi();
        Vy = zfinderM.boson().rapidity();
        Veta = zfinderM.boson().eta();
	    }

	    Jets jets = apply<JetAlg>(event, "jets").jetsByPt(Cuts::pT>20*GeV && Cuts::absrap < 2.5);
      idiscardIfAnyDeltaRLess(jets, els, 0.4);
      idiscardIfAnyDeltaRLess(jets, mus, 0.4);

	    Jets btagged;
	    const Particles allBs = apply<HeavyHadrons>(event, "HFHadrons").bHadrons(5.0*GeV);
	    Particles matchedBs;

	    for (const Jet& j : jets) {
        Jet closest_j;
        Particle closest_b;
        double minDR_j_b = 10;

        for (const Particle& bHad : allBs) {
          bool alreadyMatched = false;
          for (const Particle& bMatched : matchedBs) {
            alreadyMatched |= bMatched.isSame(bHad);
          }
          if(alreadyMatched)  continue;

          double DR_j_b = deltaR(j, bHad);
          if ( DR_j_b <= 0.3 && DR_j_b < minDR_j_b) {
            minDR_j_b = DR_j_b;
            closest_j = j;
            closest_b = bHad;
          }
        }

        if(minDR_j_b < 0.3) {
          btagged += closest_j;
          matchedBs += closest_b;
        }
	    }
	    //size_t njets = jets.size();
	    size_t ntags = btagged.size();
	    if (ntags < 1) vetoEvent;

	    _h["ib_nBJets"]->fill(1); //inclusive 1-b

	    double dYVb = fabs(Vy - btagged[0].rap());
	    double dEtaVb = fabs(Veta - btagged[0].eta());
	    double dPhiVb = deltaPhi(Vphi, btagged[0]);
	    double dRVb = sqrt(dEtaVb*dEtaVb + dPhiVb*dPhiVb);

       _h["i1b_ZpT"]   ->fill(Vpt/GeV);
       _h["i1b_ZY"]    ->fill(fabs(Vy));
       _h["i1b_dPhiZb"]->fill(dPhiVb);
       _h["i1b_dRZb"]->fill(dRVb);
       _h["i1b_dYZb"]->fill(dYVb);
       _h["i1b_bpT"]->fill(btagged[0].pt()/GeV);
       _h["i1b_bY"]->fill(btagged[0].absrap());

	    if ( ntags>1 ) {
        _h["ib_nBJets"]->fill(2); //inclusive 2-b

        double dYbb   = fabs(btagged[0].rap() - btagged[1].rap());
        double dPhibb = deltaPhi(btagged[0], btagged[1]);
        double dRbb   = deltaR(btagged[0], btagged[1]);
        double Mbb    = (btagged[0].mom() + btagged[1].mom()).mass()/GeV;
        double Ptbb   = (btagged[0].mom() + btagged[1].mom()).pt()/GeV;

        _h["i2b_ZpT"]->fill(Vpt);
        _h["i2b_dPhibb"]->fill(dPhibb);
        _h["i2b_dRbb"]->fill(dRbb);
        _h["i2b_dYbb"]->fill(dYbb);
        _h["i2b_Mbb"]->fill(Mbb);
        _h["i2b_pTbb"]->fill(Ptbb);
        _h["i2b_pTOnMbb"]->fill(Ptbb/Mbb);

      }
    }


    void finalize() {
      // routine accepts both Z->ee and Z->mm
      // data corresponds to average
      const double sf = _mode? 1.0 : 0.5;
      scale(_h, sf * crossSectionPerEvent());
    }

    protected:

      size_t _mode;

    private:

      map<string, Histo1DPtr> _h;

  };

  DECLARE_RIVET_PLUGIN(ATLAS_2020_I1788444);
}

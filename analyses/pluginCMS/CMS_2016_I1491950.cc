// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class CMS_2016_I1491950 : public Analysis {
  public:

    /// Constructor
    CMS_2016_I1491950()
      : Analysis("CMS_2016_I1491950")
    {    }

    /// Book histograms and initialise projections before the run
    void init() {

    FinalState fs(Cuts::pT > 0. && Cuts::abseta < 6.);
    PromptFinalState prompt_fs(fs);
    prompt_fs.acceptMuonDecays(true);
    prompt_fs.acceptTauDecays(true);

    // Projection for dressed electrons and muons
    Cut leptonCuts = Cuts::abseta < 2.5 and Cuts::pt > 30.*GeV;
    SpecialDressedLeptons dressedleptons(prompt_fs, leptonCuts);
    declare(dressedleptons, "DressedLeptons");

    // Neutrinos
    IdentifiedFinalState neutrinos(prompt_fs);
    neutrinos.acceptNeutrinos();
    declare(neutrinos, "Neutrinos");

    // Projection for jets
    VetoedFinalState fsForJets(fs);
    fsForJets.addVetoOnThisFinalState(dressedleptons);
    fsForJets.addVetoOnThisFinalState(neutrinos);
    declare(FastJets(fsForJets, FastJets::ANTIKT, 0.4, JetAlg::Muons::DECAY, JetAlg::Invisibles::DECAY), "Jets");

    //book hists
    book(_hist_thadpt, "d01-x02-y01");
    book(_hist_thady, "d03-x02-y01");
    book(_hist_tleppt, "d05-x02-y01");
    book(_hist_tlepy, "d07-x02-y01");
    book(_hist_ttpt, "d09-x02-y01");
    book(_hist_tty, "d13-x02-y01");
    book(_hist_ttm, "d11-x02-y01");
    book(_hist_njet, "d15-x02-y01");
    book(_hist_njets_thadpt_1, "d17-x02-y01");
    book(_hist_njets_thadpt_2, "d18-x02-y01");
    book(_hist_njets_thadpt_3, "d19-x02-y01");
    book(_hist_njets_thadpt_4, "d20-x02-y01");
    book(_hist_njets_ttpt_1, "d22-x02-y01");
    book(_hist_njets_ttpt_2, "d23-x02-y01");
    book(_hist_njets_ttpt_3, "d24-x02-y01");
    book(_hist_njets_ttpt_4, "d25-x02-y01");
    book(_hist_thady_thadpt_1, "d27-x02-y01");
    book(_hist_thady_thadpt_2, "d28-x02-y01");
    book(_hist_thady_thadpt_3, "d29-x02-y01");
    book(_hist_thady_thadpt_4, "d30-x02-y01");
    book(_hist_ttm_tty_1, "d32-x02-y01");
    book(_hist_ttm_tty_2, "d33-x02-y01");
    book(_hist_ttm_tty_3, "d34-x02-y01");
    book(_hist_ttm_tty_4, "d35-x02-y01");
    book(_hist_ttpt_ttm_1, "d37-x02-y01");
    book(_hist_ttpt_ttm_2, "d38-x02-y01");
    book(_hist_ttpt_ttm_3, "d39-x02-y01");
    book(_hist_ttpt_ttm_4, "d40-x02-y01");
    book(_histnorm_thadpt, "d42-x02-y01");
    book(_histnorm_thady, "d44-x02-y01");
    book(_histnorm_tleppt, "d46-x02-y01");
    book(_histnorm_tlepy, "d48-x02-y01");
    book(_histnorm_ttpt, "d50-x02-y01");
    book(_histnorm_tty, "d54-x02-y01");
    book(_histnorm_ttm, "d52-x02-y01");
    book(_histnorm_njet, "d56-x02-y01");
    book(_histnorm_njets_thadpt_1, "d58-x02-y01");
    book(_histnorm_njets_thadpt_2, "d59-x02-y01");
    book(_histnorm_njets_thadpt_3, "d60-x02-y01");
    book(_histnorm_njets_thadpt_4, "d61-x02-y01");
    book(_histnorm_njets_ttpt_1, "d63-x02-y01");
    book(_histnorm_njets_ttpt_2, "d64-x02-y01");
    book(_histnorm_njets_ttpt_3, "d65-x02-y01");
    book(_histnorm_njets_ttpt_4, "d66-x02-y01");
    book(_histnorm_thady_thadpt_1, "d68-x02-y01");
    book(_histnorm_thady_thadpt_2, "d69-x02-y01");
    book(_histnorm_thady_thadpt_3, "d70-x02-y01");
    book(_histnorm_thady_thadpt_4, "d71-x02-y01");
    book(_histnorm_ttm_tty_1, "d73-x02-y01");
    book(_histnorm_ttm_tty_2, "d74-x02-y01");
    book(_histnorm_ttm_tty_3, "d75-x02-y01");
    book(_histnorm_ttm_tty_4, "d76-x02-y01");
    book(_histnorm_ttpt_ttm_1, "d78-x02-y01");
    book(_histnorm_ttpt_ttm_2, "d79-x02-y01");
    book(_histnorm_ttpt_ttm_3, "d80-x02-y01");
    book(_histnorm_ttpt_ttm_4, "d81-x02-y01");

   }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // leptons
      const SpecialDressedLeptons& dressedleptons_proj = applyProjection<SpecialDressedLeptons>(event, "DressedLeptons");
      std::vector<DressedLepton> dressedLeptons = dressedleptons_proj.dressedLeptons();
      if(dressedLeptons.size() != 1) return;

      // neutrinos
      const Particles neutrinos = apply<FinalState>(event, "Neutrinos").particlesByPt();
      _nusum = FourMomentum(0., 0., 0., 0.);
      for(const Particle& neutrino : neutrinos)
          _nusum += neutrino.momentum();
      _wl = _nusum + dressedLeptons[0].momentum();

      // jets
      Cut jet_cut = (Cuts::abseta < 2.5) and (Cuts::pT > 25.*GeV);
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(jet_cut);
      Jets allJets;
      for (const Jet& jet : jets) {
        allJets.push_back(jet);
      }
      Jets bJets;
      for (const Jet& jet : allJets) {
        if (jet.bTagged()) bJets.push_back(jet);
      }

      if(bJets.size() < 2 || allJets.size() < 4) return;

      //construct top quark proxies
      double Kmin = numeric_limits<double>::max();
      for(const Jet& itaj : allJets) {
          for(const Jet& itbj : allJets) {
              if (itaj.momentum() == itbj.momentum()) continue;
              FourMomentum wh(itaj.momentum() + itbj.momentum());
              for(const Jet& ithbj : bJets) {
                  if(itaj.momentum() == ithbj.momentum() || itbj.momentum() == ithbj.momentum()) continue;
                  FourMomentum th(wh + ithbj.momentum());
                  for(const Jet& itlbj : bJets) {
                      if(itaj.momentum() == itlbj.momentum() || itbj.momentum() == itlbj.momentum() || ithbj.momentum() == itlbj.momentum()) continue;
                      FourMomentum tl(_wl + itlbj.momentum());

                      double K = pow(wh.mass() - 80.4, 2) + pow(th.mass() - 172.5, 2) + pow(tl.mass() - 172.5, 2);
                      if(K < Kmin)
                        {
                          Kmin = K;
                          _tl = tl;
                          _th = th;
                          _wh = wh;
                        }
                    }
                }
            }
        }

    _hist_thadpt->fill(_th.pt());
    _hist_thady->fill(abs(_th.rapidity()) );
    _hist_tleppt->fill(_tl.pt() );
    _hist_tlepy->fill(abs(_tl.rapidity()) );
    _histnorm_thadpt->fill(_th.pt());
    _histnorm_thady->fill(abs(_th.rapidity()) );
    _histnorm_tleppt->fill(_tl.pt() );
    _histnorm_tlepy->fill(abs(_tl.rapidity()) );
    FourMomentum tt(_tl+_th);
    _hist_ttpt->fill(tt.pt() );
    _hist_tty->fill(abs(tt.rapidity()) );
    _hist_ttm->fill(tt.mass() );
    _hist_njet->fill(min(allJets.size()-4., 4.));
    _histnorm_ttpt->fill(tt.pt() );
    _histnorm_tty->fill(abs(tt.rapidity()) );
    _histnorm_ttm->fill(tt.mass() );
    _histnorm_njet->fill(min(allJets.size()-4., 4.));
    if(allJets.size() == 4)
    {
     _hist_njets_thadpt_1->fill(_th.pt());
     _hist_njets_ttpt_1->fill(tt.pt());
     _histnorm_njets_thadpt_1->fill(_th.pt());
     _histnorm_njets_ttpt_1->fill(tt.pt());
    }
    else if(allJets.size() == 5)
    {
     _hist_njets_thadpt_2->fill(_th.pt());
     _hist_njets_ttpt_2->fill(tt.pt());
     _histnorm_njets_thadpt_2->fill(_th.pt());
     _histnorm_njets_ttpt_2->fill(tt.pt());
    }
    else if(allJets.size() == 6)
    {
     _hist_njets_thadpt_3->fill(_th.pt());
     _hist_njets_ttpt_3->fill(tt.pt());
     _histnorm_njets_thadpt_3->fill(_th.pt());
     _histnorm_njets_ttpt_3->fill(tt.pt());
    }
    else //>= 4 jets
    {
     _hist_njets_thadpt_4->fill(_th.pt());
     _hist_njets_ttpt_4->fill(tt.pt());
     _histnorm_njets_thadpt_4->fill(_th.pt());
     _histnorm_njets_ttpt_4->fill(tt.pt());
    }

    if(abs(_th.rapidity()) < 0.5)
    {
     _hist_thady_thadpt_1->fill(_th.pt());
     _histnorm_thady_thadpt_1->fill(_th.pt());
    }
    else if(abs(_th.rapidity()) < 1.0)
    {
     _hist_thady_thadpt_2->fill(_th.pt());
     _histnorm_thady_thadpt_2->fill(_th.pt());
    }
    else if(abs(_th.rapidity()) < 1.5)
    {
     _hist_thady_thadpt_3->fill(_th.pt());
     _histnorm_thady_thadpt_3->fill(_th.pt());
    }
    else if(abs(_th.rapidity()) < 2.5)
    {
     _hist_thady_thadpt_4->fill(_th.pt());
     _histnorm_thady_thadpt_4->fill(_th.pt());
    }

      if(tt.mass() >= 300. && tt.mass() < 450.)
        {
          _hist_ttm_tty_1->fill(abs(tt.rapidity()));
          _histnorm_ttm_tty_1->fill(abs(tt.rapidity()));
        }
      else if(tt.mass() >= 450. && tt.mass() < 625.)
        {
          _hist_ttm_tty_2->fill(abs(tt.rapidity()));
          _histnorm_ttm_tty_2->fill(abs(tt.rapidity()));
        }
      else if(tt.mass() >= 625. && tt.mass() < 850.)
        {
          _hist_ttm_tty_3->fill(abs(tt.rapidity()));
          _histnorm_ttm_tty_3->fill(abs(tt.rapidity()));
        }
      else if(tt.mass() >= 850. && tt.mass() < 2000.)
        {
          _hist_ttm_tty_4->fill(abs(tt.rapidity()));
          _histnorm_ttm_tty_4->fill(abs(tt.rapidity()));
        }

      if(tt.pt() < 35.)
        {
          _hist_ttpt_ttm_1->fill(tt.mass());
          _histnorm_ttpt_ttm_1->fill(tt.mass());
        }
      else if(tt.pt() < 80.)
        {
          _hist_ttpt_ttm_2->fill(tt.mass());
          _histnorm_ttpt_ttm_2->fill(tt.mass());
        }
      else if(tt.pt() < 140.)
        {
          _hist_ttpt_ttm_3->fill(tt.mass());
          _histnorm_ttpt_ttm_3->fill(tt.mass());
        }
      else if(tt.pt() < 500.)
        {
          _hist_ttpt_ttm_4->fill(tt.mass());
          _histnorm_ttpt_ttm_4->fill(tt.mass());
        }

    }


    /// Normalise histograms etc., after the run
    void finalize()
    {
      scale(_hist_thadpt, crossSection()/sumOfWeights());
      scale(_hist_thady, crossSection()/sumOfWeights());
      scale(_hist_tleppt, crossSection()/sumOfWeights());
      scale(_hist_tlepy, crossSection()/sumOfWeights());
      scale(_hist_ttpt, crossSection()/sumOfWeights());
      scale(_hist_tty, crossSection()/sumOfWeights());
      scale(_hist_ttm, crossSection()/sumOfWeights());
      scale(_hist_njet, crossSection()/sumOfWeights());
      scale(_hist_njets_thadpt_1, crossSection()/sumOfWeights());
      scale(_hist_njets_thadpt_2, crossSection()/sumOfWeights());
      scale(_hist_njets_thadpt_3, crossSection()/sumOfWeights());
      scale(_hist_njets_thadpt_4, crossSection()/sumOfWeights());
      scale(_hist_njets_ttpt_1, crossSection()/sumOfWeights());
      scale(_hist_njets_ttpt_2, crossSection()/sumOfWeights());
      scale(_hist_njets_ttpt_3, crossSection()/sumOfWeights());
      scale(_hist_njets_ttpt_4, crossSection()/sumOfWeights());
      scale(_hist_thady_thadpt_1, crossSection()/sumOfWeights()/0.5);
      scale(_hist_thady_thadpt_2, crossSection()/sumOfWeights()/0.5);
      scale(_hist_thady_thadpt_3, crossSection()/sumOfWeights()/0.5);
      scale(_hist_thady_thadpt_4, crossSection()/sumOfWeights()/1.0);
      scale(_hist_ttm_tty_1, crossSection()/sumOfWeights()/150.);
      scale(_hist_ttm_tty_2, crossSection()/sumOfWeights()/175.);
      scale(_hist_ttm_tty_3, crossSection()/sumOfWeights()/225.);
      scale(_hist_ttm_tty_4, crossSection()/sumOfWeights()/1150.);
      scale(_hist_ttpt_ttm_1, crossSection()/sumOfWeights()/35.);
      scale(_hist_ttpt_ttm_2, crossSection()/sumOfWeights()/45.);
      scale(_hist_ttpt_ttm_3, crossSection()/sumOfWeights()/60.);
      scale(_hist_ttpt_ttm_4, crossSection()/sumOfWeights()/360.);

      scale(_histnorm_thadpt, 1./_histnorm_thadpt->sumW(false));
      scale(_histnorm_thady, 1./_histnorm_thady->sumW(false));
      scale(_histnorm_tleppt, 1./_histnorm_tleppt->sumW(false));
      scale(_histnorm_tlepy, 1./_histnorm_tlepy->sumW(false));
      scale(_histnorm_ttpt, 1./_histnorm_ttpt->sumW(false));
      scale(_histnorm_tty, 1./_histnorm_tty->sumW(false));
      scale(_histnorm_ttm, 1./_histnorm_ttm->sumW(false));
      scale(_histnorm_njet, 1./_histnorm_njet->sumW(false));
      double sum_njets_thadpt = _histnorm_njets_thadpt_1->sumW(false) + _histnorm_njets_thadpt_2->sumW(false) + _histnorm_njets_thadpt_3->sumW(false) + _histnorm_njets_thadpt_4->sumW(false);
      scale(_histnorm_njets_thadpt_1, 1./sum_njets_thadpt);
      scale(_histnorm_njets_thadpt_2, 1./sum_njets_thadpt);
      scale(_histnorm_njets_thadpt_3, 1./sum_njets_thadpt);
      scale(_histnorm_njets_thadpt_4, 1./sum_njets_thadpt);
      double sum_njets_ttpt = _histnorm_njets_ttpt_1->sumW(false) + _histnorm_njets_ttpt_2->sumW(false) + _histnorm_njets_ttpt_3->sumW(false) + _histnorm_njets_ttpt_4->sumW(false);
      scale(_histnorm_njets_ttpt_1, 1./sum_njets_ttpt);
      scale(_histnorm_njets_ttpt_2, 1./sum_njets_ttpt);
      scale(_histnorm_njets_ttpt_3, 1./sum_njets_ttpt);
      scale(_histnorm_njets_ttpt_4, 1./sum_njets_ttpt);
      double sum_thady_thadpt = _histnorm_thady_thadpt_1->sumW(false) + _histnorm_thady_thadpt_2->sumW(false) + _histnorm_thady_thadpt_3->sumW(false) + _histnorm_thady_thadpt_4->sumW(false);
      scale(_histnorm_thady_thadpt_1, 1./sum_thady_thadpt/0.5);
      scale(_histnorm_thady_thadpt_2, 1./sum_thady_thadpt/0.5);
      scale(_histnorm_thady_thadpt_3, 1./sum_thady_thadpt/0.5);
      scale(_histnorm_thady_thadpt_4, 1./sum_thady_thadpt/1.0);
      double sum_ttm_tty = _histnorm_ttm_tty_1->sumW(false) + _histnorm_ttm_tty_2->sumW(false) + _histnorm_ttm_tty_3->sumW(false) + _histnorm_ttm_tty_4->sumW(false);
      scale(_histnorm_ttm_tty_1, 1./sum_ttm_tty/150.);
      scale(_histnorm_ttm_tty_2, 1./sum_ttm_tty/175.);
      scale(_histnorm_ttm_tty_3, 1./sum_ttm_tty/225.);
      scale(_histnorm_ttm_tty_4, 1./sum_ttm_tty/1150.);
      double sum_ttpt_ttm = _histnorm_ttpt_ttm_1->sumW(false) + _histnorm_ttpt_ttm_2->sumW(false) + _histnorm_ttpt_ttm_3->sumW(false) + _histnorm_ttpt_ttm_4->sumW(false);
      scale(_histnorm_ttpt_ttm_1, 1./sum_ttpt_ttm/35.);
      scale(_histnorm_ttpt_ttm_2, 1./sum_ttpt_ttm/45.);
      scale(_histnorm_ttpt_ttm_3, 1./sum_ttpt_ttm/60.);
      scale(_histnorm_ttpt_ttm_4, 1./sum_ttpt_ttm/360.);

    }


    /// @brief Special dressed lepton finder
    ///
    /// Find dressed leptons by clustering all leptons and photons
    class SpecialDressedLeptons : public FinalState {
    public:
      /// The default constructor. May specify cuts
      SpecialDressedLeptons(const FinalState& fs, const Cut& cut)
        : FinalState(cut)
      {
        setName("SpecialDressedLeptons");
        IdentifiedFinalState ifs(fs);
        ifs.acceptIdPair(PID::PHOTON);
        ifs.acceptIdPair(PID::ELECTRON);
        ifs.acceptIdPair(PID::MUON);
        declare(ifs, "IFS");
        declare(FastJets(ifs, FastJets::ANTIKT, 0.1), "LeptonJets");
      }

      /// Clone on the heap.
      virtual unique_ptr<Projection> clone() const {
        return unique_ptr<Projection>(new SpecialDressedLeptons(*this));
      }

      /// Retrieve the dressed leptons
      const vector<DressedLepton>& dressedLeptons() const { return _clusteredLeptons; }

    private:
      /// Container which stores the clustered lepton objects
      vector<DressedLepton> _clusteredLeptons;

    public:
      void project(const Event& e) {

        _theParticles.clear();
        _clusteredLeptons.clear();

        vector<DressedLepton> allClusteredLeptons;

        const Jets jets = applyProjection<FastJets>(e, "LeptonJets").jetsByPt(5.*GeV);
        for (const Jet& jet : jets) {
          Particle lepCand;
          for (const Particle& cand : jet.particles()) {
            const int absPdgId = abs(cand.pid());
            if (absPdgId == PID::ELECTRON || absPdgId == PID::MUON) {
              if (cand.pt() > lepCand.pt()) lepCand = cand;
            }
          }

          //Central lepton must be the major component
          if ((lepCand.pt() < jet.pt()/2.) || (lepCand.pid() == 0)) continue;

          DressedLepton lepton(lepCand);
          for (const Particle& cand : jet.particles()) {
            if (isSame(cand, lepCand)) continue;
            if (cand.pid() != PID::PHOTON) continue;
            lepton.addConstituent(cand, true);
          }
          allClusteredLeptons.push_back(lepton);
        }

        for (const DressedLepton& lepton : allClusteredLeptons) {
          if (accept(lepton)) {
            _clusteredLeptons.push_back(lepton);
            _theParticles.push_back(lepton.constituentLepton());
            _theParticles += lepton.constituentPhotons();
          }
        }
      }

    };


  private:

   FourMomentum _tl;
   FourMomentum _th;
   FourMomentum _wl;
   FourMomentum _wh;
   FourMomentum _nusum;

   Histo1DPtr _hist_thadpt;
   Histo1DPtr _hist_thady;
   Histo1DPtr _hist_tleppt;
   Histo1DPtr _hist_tlepy;
   Histo1DPtr _hist_ttpt;
   Histo1DPtr _hist_tty;
   Histo1DPtr _hist_ttm;
   Histo1DPtr _hist_njet;
   Histo1DPtr _hist_njets_thadpt_1;
   Histo1DPtr _hist_njets_thadpt_2;
   Histo1DPtr _hist_njets_thadpt_3;
   Histo1DPtr _hist_njets_thadpt_4;
   Histo1DPtr _hist_njets_ttpt_1;
   Histo1DPtr _hist_njets_ttpt_2;
   Histo1DPtr _hist_njets_ttpt_3;
   Histo1DPtr _hist_njets_ttpt_4;
   Histo1DPtr _hist_thady_thadpt_1;
   Histo1DPtr _hist_thady_thadpt_2;
   Histo1DPtr _hist_thady_thadpt_3;
   Histo1DPtr _hist_thady_thadpt_4;
   Histo1DPtr _hist_ttm_tty_1;
   Histo1DPtr _hist_ttm_tty_2;
   Histo1DPtr _hist_ttm_tty_3;
   Histo1DPtr _hist_ttm_tty_4;
   Histo1DPtr _hist_ttpt_ttm_1;
   Histo1DPtr _hist_ttpt_ttm_2;
   Histo1DPtr _hist_ttpt_ttm_3;
   Histo1DPtr _hist_ttpt_ttm_4;

   Histo1DPtr _histnorm_thadpt;
   Histo1DPtr _histnorm_thady;
   Histo1DPtr _histnorm_tleppt;
   Histo1DPtr _histnorm_tlepy;
   Histo1DPtr _histnorm_ttpt;
   Histo1DPtr _histnorm_tty;
   Histo1DPtr _histnorm_ttm;
   Histo1DPtr _histnorm_njet;
   Histo1DPtr _histnorm_njets_thadpt_1;
   Histo1DPtr _histnorm_njets_thadpt_2;
   Histo1DPtr _histnorm_njets_thadpt_3;
   Histo1DPtr _histnorm_njets_thadpt_4;
   Histo1DPtr _histnorm_njets_ttpt_1;
   Histo1DPtr _histnorm_njets_ttpt_2;
   Histo1DPtr _histnorm_njets_ttpt_3;
   Histo1DPtr _histnorm_njets_ttpt_4;
   Histo1DPtr _histnorm_thady_thadpt_1;
   Histo1DPtr _histnorm_thady_thadpt_2;
   Histo1DPtr _histnorm_thady_thadpt_3;
   Histo1DPtr _histnorm_thady_thadpt_4;
   Histo1DPtr _histnorm_ttm_tty_1;
   Histo1DPtr _histnorm_ttm_tty_2;
   Histo1DPtr _histnorm_ttm_tty_3;
   Histo1DPtr _histnorm_ttm_tty_4;
   Histo1DPtr _histnorm_ttpt_ttm_1;
   Histo1DPtr _histnorm_ttpt_ttm_2;
   Histo1DPtr _histnorm_ttpt_ttm_3;
   Histo1DPtr _histnorm_ttpt_ttm_4;   
   
 };



 // The hook for the plugin system
 DECLARE_RIVET_PLUGIN(CMS_2016_I1491950);

}


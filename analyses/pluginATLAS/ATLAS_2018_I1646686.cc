// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  /// @brief All-hadronic ttbar at 13 TeV
  class ATLAS_2018_I1646686 : public Analysis {
  public:

      /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2018_I1646686);

      /// Book histograms and initialise projections before the run
      void init() {

        // Get options particle-level only.
        _mode = 0;
        if ( getOption("TMODE") == "PARTICLE" ) _mode = 0;
        if ( getOption("TMODE") == "BOTH" ) _mode = 1;

        //histogram booking
        book(_h["inclusive"],1,1,1);
        bookHistograms("t_pt", 0, true);
        bookHistograms("t_y",  1, true);
        bookHistograms("t1_pt",         2);
        bookHistograms("t1_y",          3);
        bookHistograms("t2_pt",         4);
        bookHistograms("t2_y",          5);
        bookHistograms("tt_m",          6);
        bookHistograms("tt_pt",         7);
        bookHistograms("tt_y",          8);
        bookHistograms("tt_chi",        9);
        bookHistograms("tt_yboost",    10);
        bookHistograms("tt_pout",      11);
        bookHistograms("tt_dPhi",      12);
        bookHistograms("tt_Ht",        13);
        bookHistograms("tt_cosThStar", 14);

        // Projections
        Cut dressed_lep = (Cuts::abseta < 2.5) && (Cuts::pT >= 25*GeV);
        Cut eta_full = (Cuts::abseta < 5.0);

        // All final state particles
        FinalState fs(eta_full);

        // Get photons to dress leptons
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);

        // Projection to find the electrons
        PromptFinalState electrons(Cuts::abspid == PID::ELECTRON, true);
        DressedLeptons dressedelectrons(photons, electrons, 0.1, dressed_lep);
        declare(dressedelectrons, "elecs");
        DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full);

        // Projection to find the muons
        PromptFinalState muons(Cuts::abspid == PID::MUON, true);
        DressedLeptons dressedmuons(photons, muons, 0.1, dressed_lep);
        declare(dressedmuons, "muons");
        DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full);

        // Jet clustering.
        VetoedFinalState vfs;
        vfs.addVetoOnThisFinalState(ewdressedelectrons);
        vfs.addVetoOnThisFinalState(ewdressedmuons);

        FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::DECAY, JetAlg::Invisibles::DECAY);
        declare(jets, "jets");

        FastJets ljets(fs, FastJets::ANTIKT, 1.0, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(ljets, "ljets" );

        if (_mode != 0 ){
          PartonicTops partonTops;
          declare(partonTops, "partonicTops");
        }
      }


      void analyze(const Event& event) {

	if (_mode != 0){

          // Parton-level top quarks
          const Particles partonicTops = apply<PartonicTops>( event, "partonicTops").particlesByPt();
          FourMomentum top, tbar;
          bool foundT = false, foundTBar = false;
          for (const Particle& ptop : partonicTops) {
            const int pid = ptop.pid();
            if (pid == PID::TQUARK) {
              top = ptop.momentum();
              foundT = true;
            } else if (pid == -PID::TQUARK) {
              tbar = ptop.momentum();
              foundTBar = true;
            }
          }

	  FourMomentum t1_parton, t2_parton, ttbar_parton;
	  if ( foundT && foundTBar ) {
	    t1_parton = top.pT2() > tbar.pT2() ? top : tbar;
	    t2_parton = top.pT2() > tbar.pT2() ? tbar : top;
	    ttbar_parton = t1_parton + t2_parton;
	    
	    if ( t1_parton.pT() > 500*GeV && t2_parton.pT() > 350*GeV) {
	      
	      const double chi_parton = calcChi(t1_parton, t2_parton);
	      const double cosThetaStar_parton = abs(calcCosThetaStar(t1_parton, t2_parton));
	      if (cosThetaStar_parton == -99) {
		MSG_DEBUG("ttbar going faster than light! Vetoing event. Try turning of partonic tops?");
		vetoEvent;
	      }
	      const double pout_parton = abs(calcPout(t1_parton, t2_parton));
	      const double dPhi_parton = deltaPhi(t1_parton, t2_parton);
	      
	      const int randomChoice = rand() % 2;
	      const FourMomentum& randomTopParton = (randomChoice == 0) ? t1_parton : t2_parton;
	      
	      fillParton("t_pt", randomTopParton.pT()/GeV);
	      fillParton("t_y",  randomTopParton.absrap());
	      
	      fillParton("t1_pt", t1_parton.pT()/GeV);
	      fillParton("t1_y",  t1_parton.absrap());
	      fillParton("t2_pt", t2_parton.pT()/GeV);
	      fillParton("t2_y",  t2_parton.absrap());
	      
	      fillParton("tt_m",  ttbar_parton.mass()/TeV);
	      fillParton("tt_pt", ttbar_parton.pT()/GeV);
	      fillParton("tt_Ht", (t1_parton.pT() + t2_parton.pT())/GeV);
	      fillParton("tt_y",  ttbar_parton.absrap());

	      fillParton("tt_yboost", 0.5 * abs(t1_parton.rapidity() + t2_parton.rapidity()));
	      fillParton("tt_chi", chi_parton);
	      fillParton("tt_cosThStar", cosThetaStar_parton);
	      fillParton("tt_pout", pout_parton/GeV);
	      fillParton("tt_dPhi", dPhi_parton);
	    }
	  }
        }

        // Get and veto on dressed leptons
        const vector<DressedLepton> dressedElectrons = apply<DressedLeptons>(event, "elecs").dressedLeptons();
        const vector<DressedLepton> dressedMuons     = apply<DressedLeptons>(event, "muons").dressedLeptons();
        if (!dressedElectrons.empty()) vetoEvent;
        if (!dressedMuons.empty()) vetoEvent;

        // Get jets
        const Jets& all_jets  = apply<FastJets>( event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
        const FastJets& ljets_fj = apply<FastJets>( event, "ljets");
        const Jets all_ljets = ljets_fj.jetsByPt();

        // Trim the large-R jets
        Jets trimmedJets;
        fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));
        for (const Jet& jet : all_ljets)
          trimmedJets += ljets_fj.trimJet(jet, trimmer);
        trimmedJets = sortByPt(trimmedJets);

        // Check large-R jets
        Jets ljets;
        vector<bool> b_tagged;
        for (const Jet& jet : trimmedJets) {

          if (jet.pT() < 250 * GeV)  continue;
          if (jet.pT() > 3000 * GeV) continue;
          if (jet.mass() > jet.pT()) continue;
          if (jet.abseta() > 2.0 )   continue;

          ljets += jet;
          b_tagged += jet.bTagged();
        }

        if (all_jets.size() < 2)  vetoEvent;
        if (ljets.size() < 2)     vetoEvent;

        // Identify top and anti top, compute some event variables
        const FourMomentum ttbar = ljets[0].momentum() + ljets[1].momentum();
        const FourMomentum t1 = ljets[0].momentum();
        const FourMomentum t2 = ljets[1].momentum();

        const double chi = calcChi(t1, t2);
        const double cosThetaStar = abs(calcCosThetaStar(t1, t2));
	if (cosThetaStar == -99) {
	  MSG_DEBUG("real ttbar going faster than light! This should not happen. Vetoing event.");
	  vetoEvent;
	}
        const double pout = abs(calcPout(t1, t2));
        const double dPhi = deltaPhi(t1, t2);

        if ( t2.pT() < 350*GeV)  vetoEvent;
        if ( t1.pT() < 500*GeV)  vetoEvent;

        // b-tagging for particle done on large-R jets
        if (!(b_tagged[0] && b_tagged[1]))  vetoEvent;

        // Continues with signal region cuts
        if ( abs(t1.mass() - 172.5 * GeV) > 50*GeV )  vetoEvent;
        if ( abs(t2.mass() - 172.5 * GeV) > 50*GeV )  vetoEvent;

        _h["inclusive"]->fill(0);

        fillHistograms("t1_pt", t1.pT()/GeV);
        fillHistograms("t1_y",  t1.absrap());
        fillHistograms("t2_pt", t2.pT()/GeV);
        fillHistograms("t2_y",  t2.absrap());

        fillHistograms("tt_m",  ttbar.mass()/TeV);
        fillHistograms("tt_pt", ttbar.pT()/GeV);
        fillHistograms("tt_Ht", (t1.pT() + t2.pT())/GeV);
        fillHistograms("tt_y",  ttbar.absrap());

        fillHistograms("tt_yboost", 0.5 * abs(t1.rapidity() + t2.rapidity()));
        fillHistograms("tt_chi", chi);
        fillHistograms("tt_cosThStar", cosThetaStar);
        fillHistograms("tt_pout", pout/GeV);
        fillHistograms("tt_dPhi", dPhi);

      }


      void finalize() {
        // Normalize histograms
        const double sf = crossSection() / sumOfWeights();
        for (auto &hist : _h) {
          scale(hist.second, sf);
          if ((hist.first.find("_norm") != string::npos) && hist.second->integral(false)>0) hist.second->normalize(1.0, false);
        }
      }


      double calcChi(const FourMomentum& t1, const FourMomentum& t2) {
        double ystar = 0.5 * (t1.rapidity()-t2.rapidity());
        double chi = exp( 2 * abs(ystar));
        return chi;
      }

      double calcCosThetaStar(const FourMomentum& t1, const FourMomentum& t2) {
        FourMomentum ttbar = t1 + t2;
        LorentzTransform centreOfMassTrans;
        ttbar.setX(0);
        ttbar.setY(0);
        if (ttbar.betaVec().mod2() > 1) return -99;
        centreOfMassTrans.setBetaVec( -ttbar.betaVec() );
        FourMomentum t1_star = centreOfMassTrans.transform(t1);
        double cosThetaStar;
        if (t1_star.p3().mod2() >= 0){
          cosThetaStar = t1_star.pz()/t1_star.p3().mod();
        } else {
          return -99;
        }
        return cosThetaStar;
      }

      double calcPout(const FourMomentum& t1, const FourMomentum& t2) {
        Vector3 t1V = t1.p3();
        Vector3 t2V = t2.p3();
        Vector3 zUnit = Vector3(0., 0., 1.);
        Vector3 vPerp = zUnit.cross(t1V);

        double pout = vPerp.dot(t2V)/vPerp.mod();
        return pout;
      }


    protected:

      size_t _mode;

    private:

      map<string, Histo1DPtr> _h;

      //some functions for booking, filling and scaling the histograms
      void fillHistograms(std::string name, double value) {
        _h[name]->fill(value);
        _h[name + "_norm"]->fill(value);
      }

      void fillParton(std::string name, double value) {
        _h[name + "_parton"]->fill(value);
        _h[name + "_parton_norm"]->fill(value);
      }

      void bookHistograms(const std::string name, unsigned int index, bool onlyParton = false) {
        if (!onlyParton) {
          book(_h[name], index, 1, 1 );
          book(_h[name + "_norm"], index + 13, 1, 1 );
        }
        if (_mode != 0) {
          book(_h[name + "_parton"], index + 82, 1, 1 );
          book(_h[name + "_parton_norm"], index + 97, 1, 1 );
        }
      }

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2018_I1646686);

}

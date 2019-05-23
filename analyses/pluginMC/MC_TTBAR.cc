#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {


  class MC_TTBAR : public Analysis {
  public:

    /// Minimal constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_TTBAR);


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {

      _mode = 1; string pre = "onelep_"; // default is single-lepton decay mode
      if ( getOption("TTMODE") == "ALLHAD" ) { _mode = 0; pre = "allhad_"; }
      if ( getOption("TTMODE") == "ONELEP" ) { _mode = 1; pre = "onelep_"; }
      if ( getOption("TTMODE") == "TWOLEP" ) { _mode = 2; pre = "twolep_"; }
      if ( getOption("TTMODE") == "ANYLEP" ) { _mode = 3; pre = "anylep_"; }

      // A FinalState is used to select particles within |eta| < 4.2 and with pT
      // > 30 GeV, out of which the ChargedLeptons projection picks only the
      // electrons and muons, to be accessed later as "LFS".
      ChargedLeptons lfs(FinalState(Cuts::abseta < 4.2 && Cuts::pT > 30*GeV));
      declare(lfs, "LFS");

      // A second FinalState is used to select all particles in |eta| < 4.2,
      // with no pT cut. This is used to construct jets and measure missing
      // transverse energy.
      VetoedFinalState fs(FinalState(Cuts::abseta < 4.2));
      fs.addVetoOnThisFinalState(lfs);
      declare(FastJets(fs, FastJets::ANTIKT, 0.6), "Jets");
      declare(MissingMomentum(fs), "MissingET");

      // Booking of histograms
      _h["njets"] = bookHisto1D(pre + "jet_mult", 11, -0.5, 10.5);
      //
      _h["jet_1_pT"] = bookHisto1D(pre + "jet_1_pT", logspace(50, 20.0, 500.0));
      _h["jet_2_pT"] = bookHisto1D(pre + "jet_2_pT", logspace(50, 20.0, 400.0));
      _h["jet_3_pT"] = bookHisto1D(pre + "jet_3_pT", logspace(50, 20.0, 300.0));
      _h["jet_4_pT"] = bookHisto1D(pre + "jet_4_pT", logspace(50, 20.0, 200.0));
      _h["jet_HT"]   = bookHisto1D(pre + "jet_HT", logspace(50, 100.0, 2000.0));
      //
      _h["bjet_1_pT"] = bookHisto1D(pre + "jetb_1_pT", logspace(50, 20.0, 400.0));
      _h["bjet_2_pT"] = bookHisto1D(pre + "jetb_2_pT", logspace(50, 20.0, 300.0));
      //
      _h["ljet_1_pT"] = bookHisto1D(pre + "jetl_1_pT", logspace(50, 20.0, 400.0));
      _h["ljet_2_pT"] = bookHisto1D(pre + "jetl_2_pT", logspace(50, 20.0, 300.0));
      //
      if (_mode != 2)  _h["tt_mass"]   = bookHisto1D(pre + "tt_mass", 200, 300.0, 700.0);
      //
      if (_mode < 2) { // these rely on a hadronic W being part of the ttbar decay
        _h["W_mass"]        = bookHisto1D(pre + "W_mass", 75, 30, 180);
        _h["t_mass"]        = bookHisto1D(pre + "t_mass", 150, 130, 430);
        _h["t_mass_W_cut"]  = bookHisto1D(pre + "t_mass_W_cut", 150, 130, 430);
        _h["jetb_1_W_dR"]   = bookHisto1D(pre + "jetb_1_W_dR", 20, 0.0, 7.0);
        _h["jetb_1_W_deta"] = bookHisto1D(pre + "jetb_1_W_deta", 20, 0.0, 7.0);
        _h["jetb_1_W_dphi"] = bookHisto1D(pre + "jetb_1_W_dphi", 20, 0.0, M_PI);
      }
      //
      _h["jetb_1_jetb_2_dR"]   = bookHisto1D(pre + "jetb_1_jetb_2_dR", 20, 0.0, 7.0);
      _h["jetb_1_jetb_2_deta"] = bookHisto1D(pre + "jetb_1_jetb_2_deta", 20, 0.0, 7.0);
      _h["jetb_1_jetb_2_dphi"] = bookHisto1D(pre + "jetb_1_jetb_2_dphi", 20, 0.0, M_PI);
      _h["jetb_1_jetl_1_dR"]   = bookHisto1D(pre + "jetb_1_jetl_1_dR", 20, 0.0, 7.0);
      _h["jetb_1_jetl_1_deta"] = bookHisto1D(pre + "jetb_1_jetl_1_deta", 20, 0.0, 7.0);
      _h["jetb_1_jetl_1_dphi"] = bookHisto1D(pre + "jetb_1_jetl_1_dphi", 20, 0.0, M_PI);
      _h["jetl_1_jetl_2_dR"]   = bookHisto1D(pre + "jetl_1_jetl_2_dR", 20, 0.0, 7.0);
      _h["jetl_1_jetl_2_deta"] = bookHisto1D(pre + "jetl_1_jetl_2_deta", 20, 0.0, 7.0);
      _h["jetl_1_jetl_2_dphi"] = bookHisto1D(pre + "jetl_1_jetl_2_dphi", 20, 0.0, M_PI);
      if (_mode > 0) { // these rely on at least one leptonic decay mode
        _h["jetb_1_l_dR"]        = bookHisto1D(pre + "jetb_1_l_dR", 20, 0.0, 7.0);
        _h["jetb_1_l_deta"]      = bookHisto1D(pre + "jetb_1_l_deta", 20, 0.0, 7.0);
        _h["jetb_1_l_dphi"]      = bookHisto1D(pre + "jetb_1_l_dphi", 20, 0.0, M_PI);
        _h["jetb_1_l_mass"]      = bookHisto1D(pre + "jetb_1_l_mass", 40, 0.0, 500.0);
        if (_mode > 1) {
          _h["jetb_1_l2_dR"]       = bookHisto1D(pre + "jetb_1_l2_dR", 20, 0.0, 7.0);
          _h["jetb_1_l2_deta"]     = bookHisto1D(pre + "jetb_1_l2_deta", 20, 0.0, 7.0);
          _h["jetb_1_l2_dphi"]     = bookHisto1D(pre + "jetb_1_l2_dphi", 20, 0.0, M_PI);
          _h["jetb_1_l2_mass"]     = bookHisto1D(pre + "jetb_1_l2_mass", 40, 0.0, 500.0);
        }
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // Use the "LFS" projection to require at least one hard charged
      // lepton. This is an experimental signature for the leptonically decaying
      // W. This helps to reduce pure QCD backgrounds.
      const ChargedLeptons& lfs = apply<ChargedLeptons>(event, "LFS");
      MSG_DEBUG("Charged lepton multiplicity = " << lfs.chargedLeptons().size());
      for (const Particle& lepton : lfs.chargedLeptons()) {
        MSG_DEBUG("Lepton pT = " << lepton.pT());
      }

      size_t nLeps = lfs.chargedLeptons().size();
      bool leptonMultiFail = _mode == 3 && nLeps == 0; // non-all-hadronic
      leptonMultiFail |= _mode == 2 && nLeps != 2; // dilepton
      leptonMultiFail |= _mode == 1 && nLeps != 1; // single lepton
      leptonMultiFail |= _mode == 0 && nLeps != 0; // all-hadronic
      if (leptonMultiFail) {
        MSG_DEBUG("Event failed lepton multiplicity cut");
        vetoEvent;
      }

      // Use a missing ET cut to bias toward events with a hard neutrino from
      // the leptonically decaying W. This helps to reduce pure QCD backgrounds.
      // not applied in all-hadronic mode
      const Vector3& met = apply<MissingMomentum>(event, "MissingET").vectorMissingPt();
      MSG_DEBUG("Vector pT = " << met.mod() << " GeV");
      if (_mode > 0 && met.mod() < 30*GeV) {
        MSG_DEBUG("Event failed missing ET cut");
        vetoEvent;
      }

      // Use the "Jets" projection to check how many jets with pT > 30 GeV there are
      // remove jets overlapping with any lepton (dR < 0.3)
      // cut on jet multiplicity depending on ttbar decay mode
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets jets = discardIfAnyDeltaRLess(jetpro.jetsByPt(30*GeV), lfs.chargedLeptons(), 0.3);

      if (     _mode == 0 && jets.size() < 6)  vetoEvent; // all-hadronic
      else if (_mode == 1 && jets.size() < 4)  vetoEvent; // single lepton
      else if (_mode == 2 && jets.size() < 2)  vetoEvent; // dilepton
      else if (_mode == 3 && nLeps == 1 && jets.size() < 4)  vetoEvent; // non-allhadronic
      else if (_mode == 3 && nLeps == 2 && jets.size() < 2)  vetoEvent;
      MSG_DEBUG("Event failed jet multiplicity cut");

      // Fill histograms for inclusive jet kinematics 
      _h["njets"]->fill(jets.size(), weight);
      if (jets.size() > 0)  _h["jet_1_pT"]->fill(jets[0].pT()/GeV, weight);
      if (jets.size() > 1)  _h["jet_2_pT"]->fill(jets[1].pT()/GeV, weight);
      if (jets.size() > 2)  _h["jet_3_pT"]->fill(jets[2].pT()/GeV, weight);
      if (jets.size() > 3)  _h["jet_4_pT"]->fill(jets[3].pT()/GeV, weight);
      double ht = 0.0;
      for (const Jet& j : jets) { ht += j.pT(); }
      _h["jet_HT"]->fill(ht/GeV, weight);

      // Sort the jets into b-jets and light jets. We expect one hard b-jet from
      // each top decay, so our 4 hardest jets should include two b-jets. The
      // Jet::bTagged() method is equivalent to perfect experimental
      // b-tagging, in a generator-independent way.
      Jets bjets, ljets;
      for (const Jet& jet : jets) {
        if (jet.bTagged())  bjets += jet;
        else                ljets += jet;
      }
      MSG_DEBUG("Number of b-jets = " << bjets.size());
      MSG_DEBUG("Number of l-jets = " << ljets.size());
      if (bjets.size() != 2) {
        MSG_DEBUG("Event failed post-lepton-isolation b-tagging cut");
        vetoEvent;
      }
      if (_mode == 0 && ljets.size() < 4)  vetoEvent;
      else if (_mode == 1 && ljets.size() < 2)  vetoEvent;
      else if (_mode == 3 && nLeps == 1 && ljets.size() < 2)  vetoEvent;

      // Plot the pTs of the identified jets.
      _h["bjet_1_pT"]->fill(bjets[0].pT(), weight);
      _h["bjet_2_pT"]->fill(bjets[1].pT(), weight);
      // need to check size to cater for dileptonic mode
      if (ljets.size() > 0)  _h["ljet_1_pT"]->fill(ljets[0].pT(), weight);
      if (ljets.size() > 1)  _h["ljet_2_pT"]->fill(ljets[1].pT(), weight);


      // Try to reconstruct ttbar pair (doesn't really work in the dileptonic mode)
      FourMomentum ttpair = bjets[0].mom() + bjets[1].mom();
      if (_mode == 0) {
        ttpair += ljets[0].mom() + ljets[1].mom() + ljets[2].mom() + ljets[3].mom();
      }
      else if (nLeps < 2) {
        ttpair += ljets[0].mom() + ljets[1].mom();
        const FourMomentum lep = lfs.chargedLeptons()[0].mom();
        double pz = findZcomponent(lep, met);
        FourMomentum neutrino(sqrt(sqr(met.x()) + sqr(met.y()) + sqr(pz)), met.x(), met.y(), pz);
        ttpair += lep + neutrino;
      }
      if (nLeps < 2)  _h["tt_mass"]->fill(ttpair.mass()/GeV, weight);

      if (_mode < 2) {
        // Construct the hadronically decaying W momentum 4-vector from pairs of
        // non-b-tagged jets. The pair which best matches the W mass is used. We start
        // with an always terrible 4-vector estimate which should always be "beaten" by
        // a real jet pair.
        FourMomentum W(10*(sqrtS()>0.?sqrtS():14000.), 0, 0, 0);
        for (size_t i = 0; i < ljets.size()-1; ++i) {
          for (size_t j = i + 1; j < ljets.size(); ++j) {
            const FourMomentum Wcand = ljets[i].momentum() + ljets[j].momentum();
            MSG_TRACE(i << "," << j << ": candidate W mass = " << Wcand.mass()/GeV
                      << " GeV, vs. incumbent candidate with " << W.mass()/GeV << " GeV");
            if (fabs(Wcand.mass() - 80.4*GeV) < fabs(W.mass() - 80.4*GeV)) {
              W = Wcand;
            }
          }
        }
        MSG_DEBUG("Candidate W mass = " << W.mass() << " GeV");

        // There are two b-jets with which this can be combined to make the
        // hadronically decaying top, one of which is correct and the other is
        // not... but we have no way to identify which is which, so we construct
        // both possible top momenta and fill the histograms with both.
        const FourMomentum t1 = W + bjets[0].momentum();
        const FourMomentum t2 = W + bjets[1].momentum();
        _h["W_mass"]->fill(W.mass(), weight);
        _h["t_mass"]->fill(t1.mass(), weight);
        _h["t_mass"]->fill(t2.mass(), weight);

        // Placing a cut on the well-known W mass helps to reduce backgrounds
        // only done for all-hadronic and semileptonic mode (since W is hadronic)
        if (!inRange(W.mass()/GeV, 75.0, 85.0))  vetoEvent;
        MSG_DEBUG("W found with mass " << W.mass()/GeV << " GeV");

        _h["t_mass_W_cut"]->fill(t1.mass(), weight);
        _h["t_mass_W_cut"]->fill(t2.mass(), weight);

        _h["jetb_1_W_dR"]->fill(deltaR(bjets[0].momentum(), W),weight);
        _h["jetb_1_W_deta"]->fill(fabs(bjets[0].eta()-W.eta()),weight);
        _h["jetb_1_W_dphi"]->fill(deltaPhi(bjets[0].momentum(),W),weight);
      }

      _h["jetb_1_jetb_2_dR"]->fill(deltaR(bjets[0].momentum(), bjets[1].momentum()),weight);
      _h["jetb_1_jetb_2_deta"]->fill(fabs(bjets[0].eta()-bjets[1].eta()),weight);
      _h["jetb_1_jetb_2_dphi"]->fill(deltaPhi(bjets[0].momentum(),bjets[1].momentum()),weight);

      if (ljets.size() > 0) {
        _h["jetb_1_jetl_1_dR"]->fill(deltaR(bjets[0].momentum(), ljets[0].momentum()),weight);
        _h["jetb_1_jetl_1_deta"]->fill(fabs(bjets[0].eta()-ljets[0].eta()),weight);
        _h["jetb_1_jetl_1_dphi"]->fill(deltaPhi(bjets[0].momentum(),ljets[0].momentum()),weight);
        if (ljets.size() > 1) {
          _h["jetl_1_jetl_2_dR"]->fill(deltaR(ljets[0].momentum(), ljets[1].momentum()),weight);
          _h["jetl_1_jetl_2_deta"]->fill(fabs(ljets[0].eta()-ljets[1].eta()),weight);
          _h["jetl_1_jetl_2_dphi"]->fill(deltaPhi(ljets[0].momentum(),ljets[1].momentum()),weight);
        }
      }

      // lepton-centric plots
      if (_mode > 0) {
        FourMomentum l=lfs.chargedLeptons()[0].momentum();
        _h["jetb_1_l_dR"]->fill(deltaR(bjets[0].momentum(), l),weight);
        _h["jetb_1_l_deta"]->fill(fabs(bjets[0].eta()-l.eta()),weight);
        _h["jetb_1_l_dphi"]->fill(deltaPhi(bjets[0].momentum(),l),weight);
        _h["jetb_1_l_mass"]->fill(FourMomentum(bjets[0].momentum()+l).mass(), weight);

        if (nLeps > 1) {
          FourMomentum l=lfs.chargedLeptons()[1].momentum();
          _h["jetb_1_l2_dR"]->fill(deltaR(bjets[0].momentum(), l),weight);
          _h["jetb_1_l2_deta"]->fill(fabs(bjets[0].eta()-l.eta()),weight);
          _h["jetb_1_l2_dphi"]->fill(deltaPhi(bjets[0].momentum(),l),weight);
          _h["jetb_1_l2_mass"]->fill(FourMomentum(bjets[0].momentum()+l).mass(), weight);
        }
      }

    }

    double findZcomponent(const FourMomentum& lepton, const Vector3& met) const {
      // estimate z-component of momentum given lepton 4-vector and MET 3-vector
      double pz_estimate;
      double m_W = 80.399*GeV;
      double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.x() + lepton.py() * met.y());
      double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
      double b = -2*k*lepton.pz();
      double c = sqr( lepton.E() ) * sqr( met.perp() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
      if (discriminant < 0)  pz_estimate = - b / (2 * a); //if the discriminant is negative
      else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
        double absquad[2];
        for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
        if (absquad[0] < absquad[1])  pz_estimate = quad[0];
        else                          pz_estimate = quad[1];
      }
      return pz_estimate;
    }

    void finalize() {
      const double sf = crossSection() / sumOfWeights();
      for (auto hist : _h) { scale(hist.second, sf); }
    }

    //@}

  protected:

      size_t _mode;


  private:

    // @name Histogram data members
    //@{
    map<string, Histo1DPtr> _h;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTBAR);
}

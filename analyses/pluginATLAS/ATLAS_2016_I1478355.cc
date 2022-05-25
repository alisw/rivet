// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Measurement of the bbar dijet dijet cross section at 7 TeV
  class ATLAS_2016_I1478355 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1478355);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs(Cuts::abseta < 3.2);
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      fj.useInvisibles();
      declare(fj, "Jets");
      declare(HeavyHadrons(Cuts::abseta < 3.2 && Cuts::pT > 5*GeV), "BHadrons");

      book(_h["m_bb"], 1, 1, 1);
      book(_h["Delta_phi"], 2, 1, 1);
      book(_h["y_diff"], 3, 1, 1);
      book(_h["Delta_R"], 4, 1, 1);
      book(_h["pT_bb"], 5, 1, 1);
      book(_h["y_B"], 6, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(20*GeV);
      const Particles& bHadrons = apply<HeavyHadrons>(event, "BHadrons").bHadrons();

      if (jets.size() < 1)  vetoEvent;
      if (jets[0].pT() > 270*GeV){

        Jets foundBjets;
        for (const Jet& j : jets){
          for (const Particle& b : bHadrons){
            if ((deltaR(j, b) < 0.3) && (j.pT() > 20*GeV) && (fabs(j.eta()) < 2.5)) foundBjets.push_back(j);
          }
        }

        Jet firstBjet;
        Jet secondBjet;
        bool bJetPair = false;
        for (const Jet& b1 : foundBjets){
          for (const Jet& b2 : foundBjets){
            if(deltaR(b1, b2) > 0.4){
              bJetPair = true;
              firstBjet = b1;
              secondBjet = b2;
              break;
            }
          }
          if (bJetPair) break;
        }

        if (bJetPair) {
          const double mass = (firstBjet.momentum() + secondBjet.momentum()).mass();
          _h["m_bb"]->fill(mass/GeV);
          const double Delta_phi = deltaPhi(firstBjet.phi(), secondBjet.phi());
          _h["Delta_phi"]->fill(Delta_phi);
          const double y_diff = fabs(firstBjet.rapidity() - secondBjet.rapidity())/2;
          _h["y_diff"]->fill(y_diff);
          const double Delta_R = deltaR(firstBjet.momentum(), secondBjet.momentum());
          _h["Delta_R"]->fill(Delta_R);
          const double pT = (firstBjet.momentum() + secondBjet.momentum()).pT();
          _h["pT_bb"]->fill(pT/GeV);
          const double y_B = (firstBjet.rapidity() + secondBjet.rapidity())/2;
          _h["y_B"]->fill(y_B);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double xsec = crossSectionPerEvent()/picobarn; 
      scale(_h, xsec);
    }

    ///@}

    /// @name Histograms
    ///@{
     map<string, Histo1DPtr> _h;
    ///@}

  };

  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1478355);
}

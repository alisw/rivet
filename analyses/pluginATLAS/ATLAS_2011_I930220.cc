// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief ATLAS inclusive b-jet pT spectrum, di-jet mass and di-jet chi
  class ATLAS_2011_I930220: public Analysis {
  public:

    ATLAS_2011_I930220()
      : Analysis("ATLAS_2011_I930220")
    {    }


    void init() {
      FinalState fs((Cuts::etaIn(-3.5, 3.5)));
      declare(fs, "FinalState");
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      fj.useInvisibles();
      declare(fj, "Jets");
      declare(HeavyHadrons(Cuts::abseta < 3.5 && Cuts::pT > 5*GeV), "BHadrons");

      double ybins[] = { 0.0, 0.3, 0.8, 1.2, 2.1 };
      for (size_t i = 0; i < 4; ++i) {
        Histo1DPtr tmp;
        _bjetpT_SV0.add(ybins[i], ybins[i+1], book(tmp, i+1, 1, 1));
      }
      book(_bjetpT_SV0_All    ,5, 1, 1);
      book(_bjetpT_pTRel      ,6, 1, 1);
      book(_dijet_mass        ,7, 1, 1);
      book(_dijet_phi         ,8, 1, 1);
      book(_dijet_chi_110_370 ,9, 1, 1);
      book(_dijet_chi_370_850 ,10, 1, 1);

      book(_chiCounter1, "_chiCounter1");
      book(_chiCounter2, "_chiCounter2");
      book(_phiCounter, "_phiCounter1");
    }


    void analyze(const Event& evt) {

      const Particles& bHadrons = apply<HeavyHadrons>(evt, "BHadrons").bHadrons();
      const Jets& jets = apply<JetAlg>(evt, "Jets").jetsByPt(15*GeV);

      FourMomentum leadingJet, subleadingJet;
      int leadJet = 0, subJet = 0;
      for (const Jet& j : jets) {
        bool hasB = false;
        for (const Particle& b : bHadrons)
          if (deltaR(j, b) < 0.3) { hasB = true; break; }

        // Identify and classify the leading and subleading jets
        if (j.absrap() < 2.1) { ///< Move this into the jets defn
          if (!leadJet) {
            leadingJet = j.momentum();
            leadJet = (hasB && j.pT() > 40*GeV) ? 2 : 1;
          }
          else if (leadJet && !subJet) {
            subleadingJet = j.momentum();
            subJet = (hasB && j.pT() > 40*GeV) ? 2 : 1;
          }
          if (hasB) {
            _bjetpT_SV0.fill(j.absrap(), j.pT()/GeV);
            _bjetpT_SV0_All->fill(j.pT()/GeV);
            _bjetpT_pTRel->fill(j.pT()/GeV);
          }
        }
      }

      // Di-b-jet plots require both the leading and subleading jets to be b-tagged and have pT > 40 GeV
      if (leadJet == 2 && subJet == 2) {
        const double mass = FourMomentum( leadingJet + subleadingJet ).mass();
        _dijet_mass->fill(mass/GeV);

        // Plot dphi for high-mass di-b-jets
        if (mass > 110*GeV) {
          _phiCounter->fill();
          const double d_phi = deltaPhi( leadingJet.phi(), subleadingJet.phi() );
          _dijet_phi->fill(fabs(d_phi));
        }

        // Plot chi for low y_boost di-b-jets (in two high-mass bins)
        const double y_boost = 0.5 * (leadingJet.rapidity() + subleadingJet.rapidity());
        const double chi = exp( fabs( leadingJet.rapidity() - subleadingJet.rapidity() ) );
        if ( fabs(y_boost) < 1.1 ) {
          if (inRange(mass/GeV, 110, 370)) {
            _chiCounter1->fill();
            _dijet_chi_110_370->fill(chi);
          } else if (inRange(mass/GeV, 370, 850)) {
            _chiCounter2->fill();
            _dijet_chi_370_850->fill(chi);
          }
        }
      }
    }


    void finalize() {
      // Normalizing to cross-section and mass
      // Additional factors represent the division by rapidity
      const double xsec = crossSectionPerEvent()/(picobarn);
      const double chiScale1 = 1 / dbl(*_chiCounter1) / 260.0;
      const double chiScale2 = 1 / dbl(*_chiCounter2) / 480.0;
      const double phiScale  = 1 / dbl(*_phiCounter);

      _bjetpT_SV0.scale(xsec/2, this);
      scale(_bjetpT_SV0_All, xsec);
      scale(_bjetpT_pTRel, xsec);
      scale(_dijet_mass, xsec);
      scale(_dijet_phi, phiScale );
      scale(_dijet_chi_110_370, chiScale1);
      scale(_dijet_chi_370_850, chiScale2);
    }


  private:

    BinnedHistogram _bjetpT_SV0;

    Histo1DPtr _bjetpT_SV0_All;
    Histo1DPtr _bjetpT_pTRel;
    Histo1DPtr _dijet_mass;
    Histo1DPtr _dijet_phi;
    Histo1DPtr _dijet_chi_110_370;
    Histo1DPtr _dijet_chi_370_850;

    CounterPtr _chiCounter1;
    CounterPtr _chiCounter2;
    CounterPtr _phiCounter;
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2011_I930220);

}

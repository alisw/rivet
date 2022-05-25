// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS 4-jet production at 7 TeV
  class CMS_2013_I1273574 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1273574()
      : Analysis("CMS_2013_I1273574")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState cnfs((Cuts::etaIn(-4.7, 4.7)));
      declare(FastJets(cnfs, FastJets::ANTIKT, 0.5), "Jets");

      // Modified to match the HEPDATA record.
      // eta of highest pT jet
      //book(_h_jetetas[0]     ,1,1,1);
      book(_h_jetetas[0]     ,6,1,1);

      // pt of the highest pT jet
      book(_h_jetpts[0]      ,2,1,1);

      //book(_h_DeltaS         ,3,1,1);
      book(_h_DeltaS         ,12,1,1);

      //book(_h_DeltaPhiSoft   ,4,1,1);
      book(_h_DeltaPhiSoft   ,10,1,1);

      //book(_h_DeltaPtRelSoft ,5,1,1);
      book(_h_DeltaPtRelSoft ,11,1,1);

      // eta and pT of 3rd highest pT jet
      //book(_h_jetetas[2]     ,6,1,1);
      //book(_h_jetpts[2]      ,7,1,1);
      book(_h_jetetas[2]     ,8,1,1);
      book(_h_jetpts[2]      ,4,1,1);

      // eta and pT of 4th highest pT jet
      //book(_h_jetetas[3]     ,8,1,1);
      //book(_h_jetpts[3]      ,9,1,1);
      book(_h_jetetas[3]     ,9,1,1);
      book(_h_jetpts[3]      ,5,1,1);

      // eta and pT of 2nd highest pT jet
      //book(_h_jetetas[1]     ,10,1,1);
      //book(_h_jetpts[1]      ,11,1,1);
      book(_h_jetetas[1]     ,7,1,1);
      book(_h_jetpts[1]      ,3,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      /// @todo Use jetsByPt(ptGtr(20*GeV) & absetaIn(4.7)), then no need for the lower loop;
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(20*GeV);
      if (jets.size() < 4) vetoEvent;

      // Ensure that there are exactly 4 jets > 20 GeV, with two above 50 GeV
      Jets hardjets, alljets;
      for (const Jet& j : jets) {
        if (j.abseta() > 4.7) continue;
        if (j.pT() > 50*GeV) hardjets.push_back(j);
        if (j.pT() > 20*GeV) alljets.push_back(j);
      }
      if (hardjets.size() < 2 || alljets.size() != 4) vetoEvent;
      const double weight = 1.0;

      // Histogram pT and eta of all 4 jets
      for (size_t i = 0; i < 4; ++i) {
        _h_jetpts[i]->fill(alljets[i].pT()/GeV, weight);
        _h_jetetas[i]->fill(alljets[i].eta(), weight);
      }

      // Create vector sums of the hard and soft pairs of jets
      const FourMomentum p12 = alljets[0].momentum() + alljets[1].momentum();
      const FourMomentum p34 = alljets[2].momentum() + alljets[3].momentum();

      // Fill the delta(phi) between the soft jets
      const double dphisoft = deltaPhi(alljets[2], alljets[3]);
      _h_DeltaPhiSoft->fill(dphisoft, weight);

      // Fill the pT balance between the soft jets
      const double ptbalanceSoft = p34.pT() / (alljets[2].pT() + alljets[3].pT());
      _h_DeltaPtRelSoft->fill(ptbalanceSoft, weight);

      // Fill the azimuthal angle difference between the two jet pairs
      const double p12p34_trans = p12.px()*p34.px() + p12.py()*p34.py();
      const double DeltaS = acos( p12p34_trans / p12.pT() / p34.pT() );
      _h_DeltaS->fill(DeltaS, weight);
    }


    /// Normalise histograms (mostly to cross-section)
    void finalize() {
      const double invlumi = crossSection()/picobarn/sumOfWeights();
      for (size_t i = 0; i < 4; ++i) {
        scale(_h_jetpts[i], invlumi);
        scale(_h_jetetas[i], invlumi);
      }
      normalize(_h_DeltaPtRelSoft);
      normalize(_h_DeltaPhiSoft);
      normalize(_h_DeltaS);
    }


  private:

    Histo1DPtr _h_jetpts[4], _h_jetetas[4];
    Histo1DPtr _h_DeltaS, _h_DeltaPhiSoft, _h_DeltaPtRelSoft;

  };


  RIVET_DECLARE_PLUGIN(CMS_2013_I1273574);

}

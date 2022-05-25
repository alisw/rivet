#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  class ATLAS_2016_I1419070 : public Analysis {
  public:

    /// Constructor
    ATLAS_2016_I1419070() : Analysis("ATLAS_2016_I1419070")
    {  }

  public:

    void init() {

      declare(FastJets(FinalState(), FastJets::ANTIKT, 0.4), "Jets");

      book(forward_500MeV ,1, 1, 1);
      book(forward_2GeV   ,2, 1, 1);
      book(forward_5GeV   ,3, 1, 1);

      book(central_500MeV ,4, 1, 1);
      book(central_2GeV   ,5, 1, 1);
      book(central_5GeV   ,6, 1, 1);

      book(diff_500MeV, "d07-x01-y01", true);
      book(diff_2GeV  , "d08-x01-y01", true);
      book(diff_5GeV  , "d09-x01-y01", true);

      book(sum_500MeV, "d10-x01-y01", true);
      book(sum_2GeV  , "d11-x01-y01", true);
      book(sum_5GeV  , "d12-x01-y01", true);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets m_goodJets = applyProjection<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.1);

      if (m_goodJets.size() < 2)        vetoEvent;
      if (m_goodJets[0].pT() < 50*GeV)  vetoEvent;
      if (m_goodJets[1].pT() < 50*GeV)  vetoEvent;
      if (fabs(1.0 - m_goodJets[0].pT()/m_goodJets[1].pT()) > 0.5)  vetoEvent;

      bool check = m_goodJets[0].abseta() < m_goodJets[1].abseta();
      int pos_f = int(check);
      int pos_c = int(!check);

      double pt500MeV_f = CalculateNCharge(m_goodJets[pos_f], 0.5);
      double pt2GeV_f   = CalculateNCharge(m_goodJets[pos_f], 2.0);
      double pt5GeV_f   = CalculateNCharge(m_goodJets[pos_f], 5.0);
      double pT_f = m_goodJets[pos_f].pT();

      double pt500MeV_c = CalculateNCharge(m_goodJets[pos_c], 0.5);
      double pt2GeV_c   = CalculateNCharge(m_goodJets[pos_c], 2.0);
      double pt5GeV_c   = CalculateNCharge(m_goodJets[pos_c], 5.0);
      double pT_c = m_goodJets[pos_c].pT();

      forward_500MeV->fill(pT_f, pt500MeV_f);
      forward_2GeV->fill(  pT_f, pt2GeV_f);
      forward_5GeV->fill(  pT_f, pt5GeV_f);

      central_500MeV->fill(pT_c, pt500MeV_c);
      central_2GeV->fill(  pT_c, pt2GeV_c);
      central_5GeV->fill(  pT_c, pt5GeV_c);
    }

    double CalculateNCharge(Jet& jet, double pTcut=0.5) {
      unsigned int ncharge = 0;
      for (const Particle& p : jet.particles()) {
        if (p.pT() < pTcut)  continue;
        if (p.charge3())  ++ncharge;
      }
      if (ncharge > 60)  ncharge = 60;
      return double(ncharge);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      if (numEvents() > 2) {
        for (unsigned int i = 0; i < forward_500MeV->numBins(); ++i) {
          ProfileBin1D bsum  = central_500MeV->bin(i) + forward_500MeV->bin(i);
          ProfileBin1D bsum2 = central_2GeV->bin(i) + forward_2GeV->bin(i);
          ProfileBin1D bsum5 = central_5GeV->bin(i) + forward_5GeV->bin(i);
          ProfileBin1D bdiff  = central_500MeV->bin(i) - forward_500MeV->bin(i);
          ProfileBin1D bdiff2 = central_2GeV->bin(i) - forward_2GeV->bin(i);
          ProfileBin1D bdiff5 = central_5GeV->bin(i) - forward_5GeV->bin(i);

          double ydiff  = central_500MeV->bin(i).effNumEntries()? central_500MeV->bin(i).mean() : 0.0;
          double ydiff2 = central_2GeV->bin(i).effNumEntries()?   central_2GeV->bin(i).mean()   : 0.0;
          double ydiff5 = central_5GeV->bin(i).effNumEntries()?   central_5GeV->bin(i).mean()   : 0.0;
          ydiff  -= forward_500MeV->bin(i).effNumEntries()? forward_500MeV->bin(i).mean() : 0.0;
          ydiff2 -= forward_2GeV->bin(i).effNumEntries()?   forward_2GeV->bin(i).mean()   : 0.0;
          ydiff5 -= forward_5GeV->bin(i).effNumEntries()?   forward_5GeV->bin(i).mean()   : 0.0;

          double yerr  = bsum.effNumEntries()  > 1.0 ?  bsum.stdErr() : 0.0;
          double yerr2 = bsum2.effNumEntries() > 1.0 ? bsum2.stdErr() : 0.0;
          double yerr5 = bsum5.effNumEntries() > 1.0 ? bsum5.stdErr() : 0.0;

	  diff_500MeV->point(i).setY(ydiff, yerr);
          diff_2GeV->point(i).setY(ydiff2, yerr2);
          diff_5GeV->point(i).setY(ydiff5, yerr5);

          sum_500MeV->point(i).setY(bsum.effNumEntries()? bsum.mean() : 0.0, yerr);
          sum_2GeV->point(i).setY(bsum2.effNumEntries()? bsum2.mean() : 0.0, yerr2);
          sum_5GeV->point(i).setY(bsum5.effNumEntries()? bsum5.mean() : 0.0, yerr5);
        }
      }

    }


  private:

    Profile1DPtr forward_500MeV;
    Profile1DPtr forward_2GeV;
    Profile1DPtr forward_5GeV;

    Profile1DPtr central_500MeV;
    Profile1DPtr central_2GeV;
    Profile1DPtr central_5GeV;

    Scatter2DPtr sum_500MeV;
    Scatter2DPtr sum_2GeV;
    Scatter2DPtr sum_5GeV;

    Scatter2DPtr diff_500MeV;
    Scatter2DPtr diff_2GeV;
    Scatter2DPtr diff_5GeV;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1419070);
}

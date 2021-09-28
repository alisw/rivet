// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"



namespace Rivet {


  class CDF_2012_NOTE10874 : public Analysis {
  public:

    CDF_2012_NOTE10874()
      : Analysis("CDF_2012_NOTE10874")
    {}


  public:

    void init() {
      const ChargedFinalState cfs((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  0.5*GeV));
      declare(cfs, "CFS");

      int isqrts = -1;
      if (fuzzyEquals(sqrtS(), 300*GeV)) isqrts = 1;
      else if (fuzzyEquals(sqrtS(), 900*GeV)) isqrts = 2;
      else if (fuzzyEquals(sqrtS(), 1960*GeV)) isqrts = 3;
      assert(isqrts >= 0);

      book(_h_nch_transverse ,1,1,isqrts);
      book(_h_ptSumDen ,2,1,isqrts);
      book(_h_avePt ,3,1,isqrts);
    }

    // Little helper function to identify Delta(phi) regions
    inline int region_index(double dphi) {
      assert(inRange(dphi, 0.0, PI, CLOSED, CLOSED));
      if (dphi < PI/3.0) return 0;
      if (dphi < 2*PI/3.0) return 1;
      return 2;
    }


    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() < 1) {
        vetoEvent;
      }

      Particles particles = cfs.particlesByPt();
      Particle p_lead = particles[0];
      const double philead = p_lead.phi();
      const double pTlead  = p_lead.pT();

      int    tNch = 0;
      double ptSum = 0.0;
      for (const Particle& p : particles) {
        const double pT = p.pT();
        const double dPhi = deltaPhi(philead, p.phi());
        const int ir = region_index(dPhi);
        if (ir==1) {
          tNch++;
          ptSum += pT;
        }
      }

      const double dEtadPhi = 4.0*PI/3.0;

      _h_nch_transverse->fill(pTlead/GeV, tNch/dEtadPhi);
      _h_ptSumDen->fill(pTlead/GeV, ptSum/dEtadPhi);

      if (tNch > 0) {
        _h_avePt->fill(pTlead/GeV, ptSum/tNch);
      }
    }


    void finalize() {
    }


  private:

    Profile1DPtr _h_nch_transverse, _h_ptSumDen, _h_avePt;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2012_NOTE10874);

}

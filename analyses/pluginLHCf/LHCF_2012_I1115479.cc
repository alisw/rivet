// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class LHCF_2012_I1115479 : public Analysis {
  public:

    LHCF_2012_I1115479()
      : Analysis("LHCF_2012_I1115479")
    {   }


  public:

    void init() {
      declare(UnstableParticles(),"UFS");

      {Histo1DPtr tmp; _binnedHistos_y_pT.add( 8.9,  9.0, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _binnedHistos_y_pT.add( 9.0,  9.2, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _binnedHistos_y_pT.add( 9.2,  9.4, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _binnedHistos_y_pT.add( 9.4,  9.6, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _binnedHistos_y_pT.add( 9.6, 10.0, book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _binnedHistos_y_pT.add(10.0, 11.0, book(tmp, 6, 1, 1));}
    }


    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      const double dphi = TWOPI;

      for (const Particle& p : ufs.particles()) {
        if (p.pid() == 111) {
          double pT = p.pT();
          double y  = p.rapidity();

          if (pT > 0.6*GeV) continue;

          const double scaled_weight = 1.0/(dphi*pT/GeV);
          _binnedHistos_y_pT.fill(y, pT/GeV, scaled_weight);
        }
      }
    }


    void finalize() {
      _binnedHistos_y_pT.scale( 1./sumOfWeights() , this);
    }

  private:

    BinnedHistogram _binnedHistos_y_pT;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(LHCF_2012_I1115479);

}

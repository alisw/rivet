// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief BABAR Lambda_c from fragmentation
  /// @author Peter Richardson
  class BABAR_2007_S6895344 : public Analysis {
  public:

    BABAR_2007_S6895344()
      : Analysis("BABAR_2007_S6895344")
    { }


    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      book(_histOff  ,1,1,1);
      book(_sigmaOff ,2,1,1);
      book(_histOn   ,3,1,1);
      book(_sigmaOn  ,4,1,1);
    }


    void analyze(const Event& e) {
      // Loop through unstable FS particles and look for charmed mesons/baryons
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");

      const Beam beamproj = apply<Beam>(e, "Beams");
      const ParticlePair& beams = beamproj.beams();
      const FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(mom_tot.betaVec());
      const double s = sqr(beamproj.sqrtS());
      const bool onresonance = fuzzyEquals(beamproj.sqrtS(), 10.58, 2E-3);

      // Particle masses from PDGlive (accessed online 16. Nov. 2009).
      for (const Particle& p : ufs.particles()) {
        // Only looking at Lambda_c
        if (p.abspid() != 4122) continue;
        MSG_DEBUG("Lambda_c found");
        const double mH2 = 5.22780; // 2.28646^2
        const double mom = FourMomentum(cms_boost.transform(p.momentum())).p();
        const double xp = mom/sqrt(s/4.0 - mH2);

        if (onresonance) {
          _histOn  ->fill(xp);
          _sigmaOn ->fill(10.58);
        } else {
          _histOff ->fill(xp);
          _sigmaOff->fill(10.54);
        }
      }
    }


    void finalize() {
      scale(_sigmaOn , 1./sumOfWeights());
      scale(_sigmaOff, 1./sumOfWeights());
      scale(_histOn  , 1./sumOfWeights());
      scale(_histOff , 1./sumOfWeights());
    }


  private:

    //@{
    // Histograms for the continuum cross sections
    Histo1DPtr _sigmaOn ;
    Histo1DPtr _sigmaOff;
    Histo1DPtr _histOn  ;
    Histo1DPtr _histOff ;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2007_S6895344);

}

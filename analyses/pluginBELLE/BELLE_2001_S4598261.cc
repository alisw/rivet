// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief BELLE pi0 spectrum at Upsilon(4S)
  /// @author Peter Richardson
  class BELLE_2001_S4598261 : public Analysis {
  public:

    BELLE_2001_S4598261()
      : Analysis("BELLE_2001_S4598261")
    { }


    void init() {
      declare(UnstableParticles(), "UFS");
      book(_histdSigDp ,1, 1, 1); // spectrum
      book(_histMult   ,2, 1, 1); // multiplicity
      book(_weightSum, "TMP/weightSum");
    }


    void analyze(const Event& e) {
      // Find the upsilons
      Particles upsilons;
      // First in unstable final state
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");
      for (const Particle& p : ufs.particles())
        if (p.pid()==300553) upsilons.push_back(p);
      // Then in whole event if fails
      if (upsilons.empty()) {
        for(ConstGenParticlePtr p: HepMCUtils::particles(e.genEvent())) {
          if (p->pdg_id() != 300553) continue;
          ConstGenVertexPtr pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            for (ConstGenParticlePtr pp: HepMCUtils::particles(pv, Relatives::PARENTS)){
              if ( p->pdg_id() == pp->pdg_id() ) {
                passed = false;
                break;
              }
            }
          }
          if (passed) upsilons.push_back(Particle(p));
        }
      }

      // Find upsilons
      for (const Particle& p : upsilons) {
        _weightSum->fill();
        // Find the neutral pions from the decay
        vector<ConstGenParticlePtr> pions;
        findDecayProducts(p.genParticle(), pions);
        const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        for (size_t ix=0; ix<pions.size(); ++ix) {
          const double pcm = cms_boost.transform(FourMomentum(pions[ix]->momentum())).p();
          _histdSigDp->fill(pcm);
        }
        _histMult->fill(0., pions.size());
      }
    }


    void finalize() {
      scale(_histdSigDp, 1./ *_weightSum);
      scale(_histMult  , 1./ *_weightSum);
    }


  private:

    //@{
    // count of weights
    CounterPtr _weightSum;
    /// Histograms
    Histo1DPtr _histdSigDp;
    Histo1DPtr _histMult;
    //@}


    void findDecayProducts(ConstGenParticlePtr p, vector<ConstGenParticlePtr>& pions) {
      ConstGenVertexPtr dv = p->end_vertex();
      for (ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
        const int id = pp->pdg_id();
        if (id == 111) {
          pions.push_back(pp);
        } else if (pp->end_vertex())
          findDecayProducts(pp, pions);
      }
    }


  };


  DECLARE_RIVET_PLUGIN(BELLE_2001_S4598261);

}

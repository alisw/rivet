// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief BABAR Xi_c baryons from fragmentation
  /// @author Peter Richardson
  class BABAR_2005_S6181155 : public Analysis {
  public:

    BABAR_2005_S6181155()
      : Analysis("BABAR_2005_S6181155")
    { }

    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      book(_histOnResonanceA, 1,1,1);
      book(_histOnResonanceB, 2,1,1);
      book(_histOffResonance, 2,1,2);
      book(_sigma,            3,1,1);
    }

    void analyze(const Event& e) {
      // Loop through unstable FS particles and look for charmed mesons/baryons
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");

      const Beam beamproj = apply<Beam>(e, "Beams");
      const ParticlePair& beams = beamproj.beams();
      const FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(mom_tot.betaVec());
      const double s = sqr(beamproj.sqrtS());

      const bool onresonance = fuzzyEquals(beamproj.sqrtS()/GeV, 10.58, 2E-3);

      for (const Particle& p : ufs.particles()) {
        // 3-momentum in CMS frame

        const double mom = cms_boost.transform(p.momentum()).vector3().mod();
        // Only looking at Xi_c^0
        if (p.abspid() != 4132 ) continue;
        MSG_DEBUG("mom = " << mom);
        // off-resonance cross section
        if (checkDecay(p.genParticle())) {
          if (onresonance) {
            _histOnResonanceA->fill(mom);
            _histOnResonanceB->fill(mom);
          }
          else {
            _histOffResonance->fill(mom,s/sqr(10.58));
            _sigma->fill(10.58);
          }
        }
      }
    }


    void finalize() {
      scale(_histOnResonanceA, crossSection()/femtobarn/sumOfWeights()*0.2);
      scale(_histOnResonanceB, crossSection()/femtobarn/sumOfWeights()*0.45);
      scale(_histOffResonance, crossSection()/femtobarn/sumOfWeights()*0.45);
      scale(_sigma           , crossSection()/femtobarn/sumOfWeights());
    }


  private:

    //@{
    /// Histograms
    Histo1DPtr _histOnResonanceA;
    Histo1DPtr _histOnResonanceB;
    Histo1DPtr _histOffResonance;
    Histo1DPtr _sigma;
    //@}

    bool checkDecay(ConstGenParticlePtr p) {
      unsigned int nstable = 0, npip = 0, npim = 0;
      unsigned int nXim = 0, nXip = 0;
      findDecayProducts(p, nstable, npip, npim, nXip, nXim);
      int id = p->pdg_id();
      // Xi_c
      if (id == 4132) {
        if (nstable == 2 && nXim == 1 && npip == 1) return true;
      }
      else if (id == -4132) {
        if (nstable == 2 && nXip == 1 && npim == 1) return true;
      }
      return false;
    }

    void findDecayProducts(ConstGenParticlePtr p,
                           unsigned int& nstable,
                           unsigned int& npip, unsigned int& npim,
                           unsigned int& nXip, unsigned int& nXim) {
      ConstGenVertexPtr dv = p->end_vertex();
      /// @todo Use better looping
      for (ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
        int id = pp->pdg_id();
        if (id==3312) {
          ++nXim;
          ++nstable;
        } else if (id == -3312) {
          ++nXip;
          ++nstable;
        } else if(id == 111 || id == 221) {
          ++nstable;
        } else if (pp->end_vertex()) {
          findDecayProducts(pp, nstable, npip, npim, nXip, nXim);
        } else {
          if     (id !=    22) ++nstable;
          if     (id ==   211) ++npip;
          else if(id ==  -211) ++npim;
        }
      }
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2005_S6181155);

}

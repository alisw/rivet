// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// J/psi production at ATLAS
  class ATLAS_2011_S9035664: public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_S9035664);


    /// @name Analysis methods
    //@{

    void init() {
      declare(UnstableParticles(), "UFS");
      book(_nonPrRapHigh    , 14, 1, 1);
      book(_nonPrRapMedHigh , 13, 1, 1);
      book(_nonPrRapMedLow  , 12, 1, 1);
      book(_nonPrRapLow     , 11, 1, 1);
      book(_PrRapHigh       , 18, 1, 1);
      book(_PrRapMedHigh    , 17, 1, 1);
      book(_PrRapMedLow     , 16, 1, 1);
      book(_PrRapLow        , 15, 1, 1);
      book(_IncRapHigh      , 20, 1, 1);
      book(_IncRapMedHigh   , 21, 1, 1);
      book(_IncRapMedLow    , 22, 1, 1);
      book(_IncRapLow       , 23, 1, 1);
    }


    void analyze(const Event& e) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");

      for (const Particle& p : ufs.particles()) {
        if (p.abspid() != 443) continue;
        ConstGenVertexPtr gv = p.genParticle()->production_vertex();
        bool nonPrompt = false;
        if (gv) {
          for (ConstGenParticlePtr pi: HepMCUtils::particles(gv, Relatives::ANCESTORS)) {
            const PdgId pid2 = pi->pdg_id();
            if (PID::isHadron(pid2) && PID::hasBottom(pid2)) {
              nonPrompt = true;
              break;
            }
          }
        }
        double absrap = p.absrap();
        double xp = p.perp();

        if (absrap<=2.4 and absrap>2.) {
          if (nonPrompt) _nonPrRapHigh->fill(xp);
          else if (!nonPrompt) _PrRapHigh->fill(xp);
          _IncRapHigh->fill(xp);
        }
        else if (absrap<=2. and absrap>1.5) {
          if (nonPrompt) _nonPrRapMedHigh->fill(xp);
          else if (!nonPrompt) _PrRapMedHigh->fill(xp);
          _IncRapMedHigh->fill(xp);
        }
        else if (absrap<=1.5 and absrap>0.75) {
          if (nonPrompt) _nonPrRapMedLow->fill(xp);
          else if (!nonPrompt) _PrRapMedLow->fill(xp);
          _IncRapMedLow->fill(xp);
        }

        else if (absrap<=0.75) {
          if (nonPrompt) _nonPrRapLow->fill(xp);
          else if (!nonPrompt) _PrRapLow->fill(xp);
          _IncRapLow->fill(xp);
        }
      }
    }


    /// Finalize
    void finalize() {
      double factor = crossSection()/nanobarn*0.0593;

      scale(_PrRapHigh      , factor/sumOfWeights());
      scale(_PrRapMedHigh   , factor/sumOfWeights());
      scale(_PrRapMedLow    , factor/sumOfWeights());
      scale(_PrRapLow       , factor/sumOfWeights());

      scale(_nonPrRapHigh   , factor/sumOfWeights());
      scale(_nonPrRapMedHigh, factor/sumOfWeights());
      scale(_nonPrRapMedLow , factor/sumOfWeights());
      scale(_nonPrRapLow    , factor/sumOfWeights());

      scale(_IncRapHigh     , 1000.*factor/sumOfWeights());
      scale(_IncRapMedHigh  , 1000.*factor/sumOfWeights());
      scale(_IncRapMedLow   , 1000.*factor/sumOfWeights());
      scale(_IncRapLow      , 1000.*factor/sumOfWeights());

    }

    //@}


  private:

    Histo1DPtr _nonPrRapHigh;
    Histo1DPtr _nonPrRapMedHigh;
    Histo1DPtr _nonPrRapMedLow;
    Histo1DPtr _nonPrRapLow;

    Histo1DPtr _PrRapHigh;
    Histo1DPtr _PrRapMedHigh;
    Histo1DPtr _PrRapMedLow;
    Histo1DPtr _PrRapLow;

    Histo1DPtr _IncRapHigh;
    Histo1DPtr _IncRapMedHigh;
    Histo1DPtr _IncRapMedLow;
    Histo1DPtr _IncRapLow;

  };


  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2011_S9035664, ATLAS_2011_I896268);

}

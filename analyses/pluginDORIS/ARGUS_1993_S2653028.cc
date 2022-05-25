// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ARGUS pi+/-, K+/- and proton/antiproton spectrum at Upsilon(4S)
  ///
  /// @author Peter Richardson
  class ARGUS_1993_S2653028 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1993_S2653028);


    void analyze(const Event& e) {
      // Find the upsilons
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==300553)) {
	_weightSum->fill();
        Particles pionsA,pionsB,protonsA,protonsB,kaons;
        // Find the decay products we want
        findDecayProducts(p, pionsA, pionsB, protonsA, protonsB, kaons);
        LorentzTransform cms_boost;
        if (p.p3().mod() > 1*MeV)
          cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        for (size_t ix = 0; ix < pionsA.size(); ++ix) {
          FourMomentum ptemp(pionsA[ix].momentum());
          FourMomentum p2 = cms_boost.transform(ptemp);
          double pcm = cms_boost.transform(ptemp).vector3().mod();
          _histPiA->fill(pcm);
        }
        _multPiA->fill(10.58,double(pionsA.size()));
        for (size_t ix = 0; ix < pionsB.size(); ++ix) {
          double pcm = cms_boost.transform(pionsB[ix].momentum()).vector3().mod();
          _histPiB->fill(pcm);
        }
        _multPiB->fill(10.58,double(pionsB.size()));
        for (size_t ix = 0; ix < protonsA.size(); ++ix) {
          double pcm = cms_boost.transform(protonsA[ix].momentum()).vector3().mod();
          _histpA->fill(pcm);
        }
        _multpA->fill(10.58,double(protonsA.size()));
        for (size_t ix = 0; ix < protonsB.size(); ++ix) {
          double pcm = cms_boost.transform(protonsB[ix].momentum()).vector3().mod();
          _histpB->fill(pcm);
        }
        _multpB->fill(10.58,double(protonsB.size()));
        for (size_t ix = 0 ;ix < kaons.size(); ++ix) {
          double pcm = cms_boost.transform(kaons[ix].momentum()).vector3().mod();
          _histKA->fill(pcm);
          _histKB->fill(pcm);
        }
        _multK->fill(10.58,double(kaons.size()));
      }
    }


    void finalize() {
      if (_weightSum->val() > 0.) {
        scale(_histPiA, 1. / *_weightSum);
        scale(_histPiB, 1. / *_weightSum);
        scale(_histKA , 1. / *_weightSum);
        scale(_histKB , 1. / *_weightSum);
        scale(_histpA , 1. / *_weightSum);
        scale(_histpB , 1. / *_weightSum);
        scale(_multPiA, 1. / *_weightSum);
        scale(_multPiB, 1. / *_weightSum);
        scale(_multK  , 1. / *_weightSum);
        scale(_multpA , 1. / *_weightSum);
        scale(_multpB , 1. / *_weightSum);
      }
    }


    void init() {
      declare(UnstableParticles(), "UFS");

      // spectra
      book(_histPiA ,1, 1, 1);
      book(_histPiB ,2, 1, 1);
      book(_histKA  ,3, 1, 1);
      book(_histKB  ,6, 1, 1);
      book(_histpA  ,4, 1, 1);
      book(_histpB  ,5, 1, 1);
      // multiplicities
      book(_multPiA , 7, 1, 1);
      book(_multPiB , 8, 1, 1);
      book(_multK   , 9, 1, 1);
      book(_multpA  ,10, 1, 1);
      book(_multpB  ,11, 1, 1);

      book(_weightSum, "TMP/weightSum");
    }


  private:

    //@{
    /// Count of weights
    CounterPtr _weightSum;
    /// Spectra
    Histo1DPtr _histPiA, _histPiB, _histKA, _histKB, _histpA, _histpB;
    /// Multiplicities
    Histo1DPtr _multPiA, _multPiB, _multK, _multpA, _multpB;
    //@}


    void findDecayProducts(Particle parent, Particles & pionsA, Particles & pionsB,
                           Particles & protonsA, Particles & protonsB, Particles & kaons) {
      int parentId = parent.pid();
      for(const Particle & p : parent.children()) {
        int id = abs(p.pid());
        if (id == PID::PIPLUS) {
          if (parentId != PID::LAMBDA && parentId != PID::K0S) {
            pionsA.push_back(p);
            pionsB.push_back(p);
          }
          else
            pionsB.push_back(p);
        }
        else if (id == PID::PROTON) {
          if (parentId != PID::LAMBDA && parentId != PID::K0S) {
            protonsA.push_back(p);
            protonsB.push_back(p);
          }
          else
            protonsB.push_back(p);
        }
        else if (id == PID::KPLUS) {
          kaons.push_back(p);
        }
        else if (!p.children().empty())
          findDecayProducts(p, pionsA, pionsB, protonsA, protonsB, kaons);
      }
    }


  };



  RIVET_DECLARE_ALIASED_PLUGIN(ARGUS_1993_S2653028, ARGUS_1993_I340894);

}

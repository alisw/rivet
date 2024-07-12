// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ARGUS vector meson production
  ///
  /// @author Peter Richardson
  class ARGUS_1993_S2789213 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1993_S2789213);


    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix) {
        for(unsigned int iy=0;iy<5;++iy) {
          std::ostringstream title;
          title << "/TMP/MULT_" << ix << "_" << iy;
          book(_mult[ix][iy],title.str());
        }
      }

      book(_hist_cont_KStarPlus , 4, 1, 1);
      book(_hist_Ups1_KStarPlus , 5, 1, 1);
      book(_hist_Ups4_KStarPlus , 6, 1, 1);

      book(_hist_cont_KStar0    , 7, 1, 1);
      book(_hist_Ups1_KStar0    , 8, 1, 1);
      book(_hist_Ups4_KStar0    , 9, 1, 1);

      book(_hist_cont_Rho0      ,10, 1, 1);
      book(_hist_Ups1_Rho0      ,11, 1, 1);
      book(_hist_Ups4_Rho0      ,12, 1, 1);

      book(_hist_cont_Omega     ,13, 1, 1);
      book(_hist_Ups1_Omega     ,14, 1, 1);


      book(_weightSum_cont,"TMP/weightSumcont");
      book(_weightSum_Ups1,"TMP/weightSumUps1");
      book(_weightSum_Ups4,"TMP/weightSumUps4");
    }


    void analyze(const Event& e) {
      // Find the upsilons
      // First in unstable final state
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 || Cuts::pid==300553);
      // continuum
      if (upsilons.empty()) {
        _weightSum_cont->fill();
        for (const Particle& p : ufs.particles()) {
          int id = p.abspid();
          double xp = 2.*p.E()/sqrtS();
          double beta = p.p3().mod()/p.E();
          if (id == 113) {
            _hist_cont_Rho0->fill(xp, 1./beta);
            _mult[0][1]->fill();
          }
          else if (id == 313) {
            _hist_cont_KStar0->fill(xp, 1./beta);
            _mult[0][2]->fill();
          }
          else if (id == 223) {
            _hist_cont_Omega->fill(xp, 1./beta);
            _mult[0][0]->fill();
          }
          else if (id == 323) {
            _hist_cont_KStarPlus->fill(xp,1./beta);
            _mult[0][3]->fill();
          }
          else if (id == 333) {
            _mult[0][4]->fill();
          }
        }
      }
      // found an upsilon
      else {
        for (const Particle& ups : upsilons) {
          const int parentId = ups.pid();
          if(parentId == 553)
            _weightSum_Ups1->fill();
          else
            _weightSum_Ups4->fill();
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups,unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 0.001)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          double mass = ups.mass();
          for( const Particle & p : unstable) {
            int id = p.abspid();
            FourMomentum p2 = cms_boost.transform(p.momentum());
            double xp = 2.*p2.E()/mass;
            double beta = p2.p3().mod()/p2.E();
            if (id == 113) {
              if (parentId == 553) {
                _hist_Ups1_Rho0->fill(xp,1./beta);
                _mult[1][1]->fill();
              }
              else {
                _hist_Ups4_Rho0->fill(xp,1./beta);
                _mult[2][1]->fill();
              }
            }
            else if (id == 313) {
              if (parentId == 553) {
                _hist_Ups1_KStar0->fill(xp,1./beta);
                _mult[1][2]->fill();
              }
              else {
                _hist_Ups4_KStar0->fill(xp,1./beta);
                _mult[2][2]->fill();
              }
            }
            else if (id == 223) {
              if (parentId == 553) {
                _hist_Ups1_Omega->fill(xp,1./beta);
                _mult[1][0]->fill();
              }
              else {
                _mult[2][0]->fill();
              }
            }
            else if (id == 323) {
              if (parentId == 553) {
                _hist_Ups1_KStarPlus->fill(xp,1./beta);
                _mult[1][3]->fill();
              }
              else {
                _hist_Ups4_KStarPlus->fill(xp,1./beta);
                _mult[2][3]->fill();
              }
            }
            else if (id == 333) {
              if (parentId == 553) {
                _mult[1][4]->fill();
              }
              else {
                _mult[2][4]->fill();
              }
            }
          }
        }
      }
    }


    void finalize() {
      // multiplicities
      vector<CounterPtr> scales = {_weightSum_cont,_weightSum_Ups1,_weightSum_Ups4};
      for(unsigned int ix=0;ix<3;++ix) {
        if(scales[ix]->val() <= 0.) continue;
        for(unsigned int iy=0;iy<5;++iy) {
          // skip Upsilon(4S) -> omega, just an upper limit
          if(ix==2&&iy==0) continue;
          Scatter2DPtr scatter;
          book(scatter,ix+1, 1, iy+1, true);
          scale(_mult[ix][iy],1./ *scales[ix]);
          scatter->point(0).setY(_mult[ix][iy]->val(),_mult[ix][iy]->err());
        }
      }
      // spectra
      if (_weightSum_cont->val() > 0.) {
        scale(_hist_cont_KStarPlus, 1. / *_weightSum_cont);
        scale(_hist_cont_KStar0   , 1. / *_weightSum_cont);
        scale(_hist_cont_Rho0     , 1. / *_weightSum_cont);
        scale(_hist_cont_Omega    , 1. / *_weightSum_cont);
      }
      if (_weightSum_Ups1->val() > 0.) {
        scale(_hist_Ups1_KStarPlus, 1. / *_weightSum_Ups1);
        scale(_hist_Ups1_KStar0   , 1. / *_weightSum_Ups1);
        scale(_hist_Ups1_Rho0     , 1. / *_weightSum_Ups1);
        scale(_hist_Ups1_Omega    , 1. / *_weightSum_Ups1);
      }
      if (_weightSum_Ups4->val() > 0.) {
        scale(_hist_Ups4_KStarPlus, 1. / *_weightSum_Ups4);
        scale(_hist_Ups4_KStar0   , 1. / *_weightSum_Ups4);
        scale(_hist_Ups4_Rho0     , 1. / *_weightSum_Ups4);
      }
    }


  private:

    Histo1DPtr _hist_cont_KStarPlus, _hist_Ups1_KStarPlus, _hist_Ups4_KStarPlus;
    Histo1DPtr _hist_cont_KStar0, _hist_Ups1_KStar0, _hist_Ups4_KStar0   ;
    Histo1DPtr _hist_cont_Rho0, _hist_Ups1_Rho0,  _hist_Ups4_Rho0;
    Histo1DPtr _hist_cont_Omega, _hist_Ups1_Omega;

    CounterPtr _mult[3][5];
    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups4;


    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
        if (id == 113 || id == 313 || id == 323 ||
            id == 333 || id == 223 ) {
          unstable.push_back(p);
        }
        else if(!p.children().empty())
          findDecayProducts(p, unstable);
      }
    }

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ARGUS_1993_S2789213, ARGUS_1993_I356616);

}

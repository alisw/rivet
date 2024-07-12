// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda in Upsilon(1S) decay and nearby continuum
  class MD1_1994_I362599 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MD1_1994_I362599);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histos
      book(_weightSum_cont,"TMP/weightSumcont");
      book(_weightSum_Ups1,"TMP/weightSumUps1");
      book(_h_spect,1,1,1);
      for (unsigned int ix=0; ix<2; ++ix) {
        for (unsigned int iy=0; iy<2; ++iy) {
          book(_mult[ix][iy],"/TMP/MULT_" +toString(ix) + "_" +toString(iy));
        }
      }
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for (const Particle & p: mother.children()) {
        const int id = p.abspid();
        if (id==PID::LAMBDA   || id==PID::XIMINUS) unstable.push_back(p);
        if (!p.children().empty()) findDecayProducts(p, unstable);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the upsilons
      // First in unstable final state
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      // continuum
      if (upsilons.empty()) {
        _weightSum_cont->fill();
        // Unstable particles
        for (const Particle& p : ufs.particles(Cuts::abspid==PID::LAMBDA)) {
          (void)p; // to suppress "unused variable" warning
          _mult[1][0]->fill();
        }
      }
      else {
        for (const Particle& ups : upsilons) {
          _weightSum_Ups1->fill();
          Particles unstable;
          LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          // Find the decay products we want
          findDecayProducts(ups,unstable);
          for(const Particle & p : unstable)  {
            int id = p.abspid();
            if(id==PID::LAMBDA) {
              double xp = 2.*boost.transform(p.momentum()).p3().mod()/ups.mass();
              _h_spect->fill(xp);
              _mult[0][0]->fill();
            }
            else if(id==PID::XIMINUS) {
              _mult[0][1]->fill();
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_weightSum_Ups1->val() > 0.) {
        scale(_h_spect,1./ *_weightSum_Ups1);
        for(unsigned int iy=0;iy<2;++iy) {
          Scatter2DPtr scatter;
          book(scatter,2+iy*2, 1, 1, true);
          if(_weightSum_Ups1->val() <= 0.) {
            scatter->point(0).setY(0.,0.);
          }
          else {
            scale(_mult[0][iy],1./ *_weightSum_Ups1);
            scatter->point(0).setY(_mult[0][iy]->val(),_mult[0][iy]->err());
          }
        }
      }
      if (_weightSum_cont->val() > 0.) {
        scale(_mult[1][0],1./ *_weightSum_cont);
        for(unsigned int ix=0;ix<2;++ix) {
          Scatter2DPtr scatter;
          book(scatter,3, 1, 1+ix, true);
          if(inRange(sqrtS(),7.2,10.) && ix==1)
            scatter->point(0).setY(_mult[1][0]->val(),_mult[1][0]->err());
          else if(inRange(sqrtS(),7.2,9.4) && ix==0)
            scatter->point(0).setY(_mult[1][0]->val(),_mult[1][0]->err());
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spect;
    CounterPtr _weightSum_cont,_weightSum_Ups1;
    CounterPtr _mult[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MD1_1994_I362599);

}

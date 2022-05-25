// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief 
  class CLEO_1985_I205668 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1985_I205668);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
      book(_weightSum_cont,"TMP/weightSumcont");
      book(_weightSum_Ups1,"TMP/weightSumUps1");
      // multiplcities
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<12;++iy) {
	  std::ostringstream title;
	  title << "/TMP/MULT_" << ix << "_" << iy;
	  book(_mult[ix][iy],title.str());
	}
      }
      // cont spectra
      book(_h_cont_pip   , 1,1,1);
      book(_h_cont_Kp    , 2,1,1);
      book(_h_cont_p     , 3,1,1);
      book(_h_cont_pi0   , 4,1,1);
      book(_h_cont_K0    , 5,1,1);
      book(_h_cont_lam   , 6,1,1);
      book(_h_cont_xi    , 7,1,1);
      book(_h_cont_rho   , 8,1,1);
      book(_h_cont_Kstarp, 9,1,1);
      book(_h_cont_Kstar0,10,1,1);
      book(_h_cont_phi   ,11,1,1);
      // ups spectra
      book(_h_ups1_pip   , 1,1,2);
      book(_h_ups1_Kp    , 2,1,2);
      book(_h_ups1_p     , 3,1,2);
      book(_h_ups1_pi0   , 4,1,2);
      book(_h_ups1_K0    , 5,1,2);
      book(_h_ups1_lam   , 6,1,2);
      book(_h_ups1_xi    , 7,1,2);
      book(_h_ups1_rho   , 8,1,2);
      book(_h_ups1_Kstarp, 9,1,2);
      book(_h_ups1_Kstar0,10,1,2);
      book(_h_ups1_phi   ,11,1,2);
    }
    
    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = p.abspid();
	if (id == PID::PIPLUS || id==PID::KPLUS   || id==PID::PROTON ||
	    id==PID::PI0      || id==PID::K0S     || id==PID::K0L    ||
	    id==PID::LAMBDA   || id==PID::XIMINUS || id==PID::RHO0   ||
	    id==323 || id==313 || id==225 || id==PID::PHI) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
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
	const FinalState & fs = apply<FinalState>(event, "FS");
	// FS particles
        for (const Particle& p : fs.particles()) {
          int id = p.abspid();
          double xp = 2.*p.p3().mod()/sqrtS();
	  if(id==PID::PIPLUS) {
	    _h_cont_pip->fill(xp);
	    _mult[1][0]->fill();
	  }
	  else if(id==PID::KPLUS) {
	    _h_cont_Kp->fill(xp);
	    _mult[1][1]->fill();
	  }
	  else if(id==PID::PROTON) {
	    _h_cont_p->fill(xp);
	    _mult[1][2]->fill();
	  }
	}
	// Unstable particles
        for (const Particle& p : ufs.particles()) {
          int id = p.abspid();
          double xp = 2.*p.p3().mod()/sqrtS();
	  if(id==PID::PI0) {
	    _h_cont_pi0->fill(xp);
	    _mult[1][3]->fill();
	  }
	  else if(id==PID::K0S || id==PID::K0L) {
	    _h_cont_K0->fill(xp);
	    _mult[1][4]->fill();
	  }
	  else if(id==PID::LAMBDA) {
	    _h_cont_lam->fill(xp);
	    _mult[1][5]->fill();
	  }
	  else if(id==PID::XIMINUS) {
	    _h_cont_xi->fill(xp);
	    _mult[1][6]->fill();
	  }
	  else if(id==PID::RHO0) {
	    _h_cont_rho->fill(xp);
	    _mult[1][7]->fill();
	  }
	  else if(id==323) {
	    _h_cont_Kstarp->fill(xp);
	    _mult[1][8]->fill();
	  }
	  else if(id==313) {
	    _h_cont_Kstar0->fill(xp);
	    _mult[1][9]->fill();
	  }
	  else if(id==PID::PHI) {
	    _h_cont_phi->fill(xp);
	    _mult[1][10]->fill();
	  }
	  else if(id==225) {
	    _mult[1][11]->fill();
	  }
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
	    double xp = 2.*boost.transform(p.momentum()).p3().mod()/ups.mass();
	    if(id==PID::PIPLUS) {
	      _h_ups1_pip->fill(xp);
	      _mult[0][0]->fill();
	    }
	    else if(id==PID::KPLUS) {
	      _h_ups1_Kp->fill(xp);
	      _mult[0][1]->fill();
	    }
	    else if(id==PID::PROTON) {
	      _h_ups1_p->fill(xp);
	      _mult[0][2]->fill();
	    }
	    else if(id==PID::PI0) {
	      _h_ups1_pi0->fill(xp);
	      _mult[0][3]->fill();
	    }
	    else if(id==PID::K0S || id==PID::K0L) {
	      _h_ups1_K0->fill(xp);
	      _mult[0][4]->fill();
	    }
	    else if(id==PID::LAMBDA) {
	      _h_ups1_lam->fill(xp);
	      _mult[0][5]->fill();
	    }
	    else if(id==PID::XIMINUS) {
	      _h_ups1_xi->fill(xp);
	      _mult[0][6]->fill();
	    }
	    else if(id==PID::RHO0) {
	      _h_ups1_rho->fill(xp);
	      _mult[0][7]->fill();
	    }
	    else if(id==323) {
	      _h_ups1_Kstarp->fill(xp);
	      _mult[0][8]->fill();
	    }
	    else if(id==313) {
	      _h_ups1_Kstar0->fill(xp);
	      _mult[0][9]->fill();
	    }
	    else if(id==PID::PHI) {
	      _h_ups1_phi->fill(xp);
	      _mult[0][10]->fill();
	    }
	    else if(id==225) {
	      _mult[0][11]->fill();
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // multiplicities
      vector<CounterPtr> scales = {_weightSum_Ups1,_weightSum_cont};
      for(unsigned int ix=0;ix<12;++ix) {
	Scatter2DPtr scatter;
	book(scatter,ix+12, 1, 1, true);
	for(unsigned int iy=0;iy<2;++iy) {
	  if(scales[iy]->val() <= 0.) {
	    scatter->point(iy).setY(0.,0.);
	  }
	  else {
	    scale(_mult[iy][ix],1./ *scales[iy]);
	    scatter->point(iy).setY(_mult[iy][ix]->val(),_mult[iy][ix]->err());
	  }
	}
      }
      // spectra
      if (_weightSum_cont->val() > 0.) {
        scale(_h_cont_pip   , 1. / *_weightSum_cont);
        scale(_h_cont_Kp    , 1. / *_weightSum_cont);
        scale(_h_cont_p     , 1. / *_weightSum_cont);
        scale(_h_cont_pi0   , 1. / *_weightSum_cont);
        scale(_h_cont_K0    , 1. / *_weightSum_cont);
        scale(_h_cont_lam   , 1. / *_weightSum_cont);
        scale(_h_cont_xi    , 1. / *_weightSum_cont);
        scale(_h_cont_rho   , 1. / *_weightSum_cont);
        scale(_h_cont_Kstarp, 1. / *_weightSum_cont);
        scale(_h_cont_Kstar0, 1. / *_weightSum_cont);
        scale(_h_cont_phi   , 1. / *_weightSum_cont);
      }
      if (_weightSum_Ups1->val() > 0.) {
        scale(_h_ups1_pip   , 1. / *_weightSum_Ups1);
        scale(_h_ups1_Kp    , 1. / *_weightSum_Ups1);
        scale(_h_ups1_p     , 1. / *_weightSum_Ups1);
        scale(_h_ups1_pi0   , 1. / *_weightSum_Ups1);
        scale(_h_ups1_K0    , 1. / *_weightSum_Ups1);
        scale(_h_ups1_lam   , 1. / *_weightSum_Ups1);
        scale(_h_ups1_xi    , 1. / *_weightSum_Ups1);
        scale(_h_ups1_rho   , 1. / *_weightSum_Ups1);
        scale(_h_ups1_Kstarp, 1. / *_weightSum_Ups1);
        scale(_h_ups1_Kstar0, 1. / *_weightSum_Ups1);
        scale(_h_ups1_phi   , 1. / *_weightSum_Ups1);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_cont_pip,_h_cont_Kp,_h_cont_p,_h_cont_pi0,_h_cont_K0,_h_cont_lam,
      _h_cont_xi,_h_cont_rho,_h_cont_Kstarp,_h_cont_Kstar0,_h_cont_phi;
    Histo1DPtr _h_ups1_pip,_h_ups1_Kp,_h_ups1_p,_h_ups1_pi0,_h_ups1_K0,_h_ups1_lam,
      _h_ups1_xi,_h_ups1_rho,_h_ups1_Kstarp,_h_ups1_Kstar0,_h_ups1_phi;
    CounterPtr _weightSum_cont,_weightSum_Ups1;
    CounterPtr _mult[2][12];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1985_I205668);

}

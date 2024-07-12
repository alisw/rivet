// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> gamma chi_c1,2
  class CLEO_2009_I832707 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2009_I832707);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::pid==20443 || Cuts::pid==445), "UFS");
      declare(FinalState(), "FS");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // types of final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p :  fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // loop over chi_c states
      Particle chi;
      bool matched = false;
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p :  ufs.particles()) {
       	if(p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
	if(ncount==1) {
	  matched = true;
	  for(auto const & val : nRes) {
	    if(val.first==PID::PHOTON) {
	      if(val.second!=1) {
	      matched = false;
	      break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    chi=p;
	    break;
	  }
	}
      }
      if(!matched) vetoEvent;
      // have chi_c find psi2S 
      if(chi.parents().empty() || chi.children().size()!=2) vetoEvent;
      Particle psi2S = chi.parents()[0];
      if(psi2S.pid()!=100443 || psi2S.children().size()!=2) vetoEvent;
      // then the first photon
      Particle gamma1;
      if(psi2S.children()[0].pid()==PID::PHOTON)
	gamma1 = psi2S.children()[0];
      else if(psi2S.children()[1].pid()==PID::PHOTON)
	gamma1 = psi2S.children()[1];
      else
	vetoEvent;
      // then the J/psi and second photon
      Particle JPsi,gamma2;
      if(chi.children()[0].pid()==PID::PHOTON &&
	 chi.children()[1].pid()==443) {
	gamma2 = chi.children()[0];
	JPsi   = chi.children()[1];
      }
      else if(chi.children()[1].pid()==PID::PHOTON &&
	      chi.children()[0].pid()==443) {
	gamma2 = chi.children()[1];
	JPsi   = chi.children()[0];
      }
      else
	vetoEvent;
      // finally the leptons from J/psi decay
      if(JPsi.children().size()!=2) vetoEvent;
      if(JPsi.children()[0].pid()!=-JPsi.children()[1].pid()) vetoEvent;
      if(JPsi.children()[0].abspid()!=PID::EMINUS &&
	 JPsi.children()[0].abspid()!=PID::MUON) vetoEvent;
      Particle lm = JPsi.children()[0];
      Particle lp = JPsi.children()[1];
      if(lm.pid()<0) swap(lm,lp);
      // type chi state
      unsigned int ichi= chi.pid()==445 ? 1 : 0;
      // axis in the chi frame
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(chi.momentum().betaVec());
      FourMomentum pGamma2 = boost1.transform(gamma2.momentum());
      Vector3 axis1 = pGamma2.p3().unit();
      // cos thetaxs distributions
      FourMomentum pJpsi = boost1.transform(JPsi.momentum());
      FourMomentum plp   = boost1.transform(  lp.momentum());
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pJpsi.betaVec());
      Vector3 axis2 = boost2.transform(plp).p3().unit();
      Vector3 e2z = gamma2.momentum().p3().unit();
      _h[ichi]->fill(abs(e2z.dot(axis2)));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2009_I832707);

}

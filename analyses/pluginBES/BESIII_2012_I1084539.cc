// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi -> eta(1405)
  class BESIII_2012_I1084539 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1084539);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      declare(Beam(), "Beams");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
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
      Particles psi = apply<UnstableParticles>(event,"UFS").particles();
      if(psi.size()!=1 || psi[0].children().size()!=2) vetoEvent;
      Particle gamma,eta;
      if(psi[0].children()[0].pid()==PID::GAMMA &&
	 psi[0].children()[1].pid()==9020221) {
	gamma = psi[0].children()[0];
	eta   = psi[0].children()[1];
      }
      else if(psi[0].children()[1].pid()==PID::GAMMA &&
	      psi[0].children()[0].pid()==9020221) {
	gamma = psi[0].children()[1];
	eta   = psi[0].children()[0];
      }
      else
	vetoEvent;
      if(eta.children().size()!=2) vetoEvent;
      // eta(1405) -> f0(980) pi0
      Particle f0,pi0;
      if(eta.children()[0].pid()==PID::PI0 &&
	 eta.children()[1].pid()==9010221) {
	pi0 = eta.children()[0];
	f0  = eta.children()[1];
      }
      else if(eta.children()[1].pid()==PID::PI0 &&
	      eta.children()[0].pid()==9010221) {
	pi0 = eta.children()[1];
	f0  = eta.children()[0];
      }
      else
	vetoEvent;
      double cTheta = axis.dot(gamma.p3().unit());
      _h[1]->fill(cTheta);
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(psi[0].momentum().betaVec());
      FourMomentum pGamma = boost1.transform(gamma.momentum());
      Vector3 axis1 = pGamma.p3().unit();
      FourMomentum pEta = boost1.transform(eta.momentum());
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pEta.betaVec());
      FourMomentum pf0 = boost2.transform(boost1.transform(f0.momentum()));
      _h[0]->fill(pf0.p3().unit().dot(axis1));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1084539);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> omega phi
  class BESIII_2022_I1900094 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I1900094);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==PID::D0), "UFS");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS || id == PID::K0S || id == PID::K0L ) {
	  ++nstable;
	}
	else if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & D0 : apply<UnstableParticles>(event, "UFS").particles()) {
	Particle omega,phi;
	if(D0.children()[0].pid()==PID::OMEGA && D0.children()[1].pid()==PID::PHI) {
	  omega = D0.children()[0];
	  phi   = D0.children()[1];
	}
	else if(D0.children()[1].pid()==PID::OMEGA && D0.children()[0].pid()==PID::PHI) {
	  omega = D0.children()[1];
	  phi   = D0.children()[0];
	}
	else
	  continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(D0.momentum().betaVec());
	// first the phi
	if(phi.children().size()==2 &&
	   phi.children()[0].abspid()==PID::KPLUS &&
	   phi.children()[1].abspid()==PID::KPLUS ) {
	  // first boost all relevant momenta to D0 rest frame
	  FourMomentum pD0  = boost.transform(D0.momentum());
	  FourMomentum pphi = boost.transform(phi.momentum());
	  FourMomentum pK   = boost.transform(phi.children()[0].momentum());
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pphi.betaVec());
	  Vector3 axis1 = boost2.transform(pD0).p3().unit();
	  Vector3 axis2 = boost2.transform(pK ).p3().unit();
	  _h[1]->fill(abs(axis1.dot(axis2)));
	}
	// then the omega
	unsigned int nstable(0);
	Particles pip, pim, pi0;
	findDecayProducts(omega, nstable, pip, pim, pi0);
	if(nstable==3 && pip.size()==1 && pip.size()==1 && pi0.size()==1) {
	  FourMomentum pD0    = boost.transform(D0.momentum());
	  FourMomentum pomega = boost.transform(omega.momentum() );
	  FourMomentum ppip   = boost.transform(pip[0].momentum());
	  FourMomentum ppim   = boost.transform(pim[0].momentum());
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pomega.betaVec());
	  Vector3 axis = boost2.transform(pD0).p3().unit();
	  Vector3 pp   = boost2.transform(ppip).p3();
	  Vector3 pm   = boost2.transform(ppim).p3();
	  Vector3 axis2 = pp.cross(pm).unit();
	  _h[0]->fill(abs(axis.dot(axis2)));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I1900094);

}

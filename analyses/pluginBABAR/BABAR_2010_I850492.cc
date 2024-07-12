// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Upsilon_2 -> pi+ pi- Upsilon
  class BABAR_2010_I850492 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2010_I850492);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==20555);
      declare(ufs, "UFS");
      DecayedParticles Upsilon2(ufs);
      Upsilon2.addStable(PID::PI0);
      Upsilon2.addStable(553);
      declare(Upsilon2, "Upsilon2");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 211,1}, {-211,1}, {553,1} };
      DecayedParticles Upsilon2 = apply<DecayedParticles>(event, "Upsilon2");
      // loop over particles
      for(unsigned int ix=0;ix<Upsilon2.decaying().size();++ix) {
	if ( !Upsilon2.modeMatches(ix,3,mode) ) continue;
       	const Particle & pip= Upsilon2.decayProducts()[ix].at( 211)[0];
       	const Particle & pim= Upsilon2.decayProducts()[ix].at(-211)[0];
       	const Particle & ups= Upsilon2.decayProducts()[ix].at( 553)[0];
	FourMomentum ptot = pip.momentum()+pim.momentum();
	_h[0]->fill(ptot.mass());
	// boost to Upsilon_2 rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(Upsilon2.decaying()[ix].momentum().betaVec());
	FourMomentum pDir = boost.transform(ptot);
	Matrix3 ptoz(-pDir.p3().unit(), Vector3(0,0,1));
	boost.preMult(ptoz);
	FourMomentum p2 = boost.transform(ups.momentum());
	FourMomentum ppip = boost.transform(pip.momentum());
	FourMomentum ppim = boost.transform(pim.momentum());
	ptot = ppip+ppim;
	// pion angle
	LorentzTransform boostPi = LorentzTransform::mkFrameTransformFromBeta(ptot.betaVec());
      	Vector3 axisPi = boostPi.transform(ppip).p3().unit();
	double cosPi = axisPi.dot(ptot.p3().unit());
      	_h[2]->fill(abs(cosPi));
	if(ups.children().size()!=2) continue;
	Particle ep,em;
	if ( ups.children()[0].pid()==-ups.children()[1].pid() &&
	     (ups.children()[0].abspid()==11 || ups.children()[0].abspid()==13)) {
	  ep = ups.children()[0];
	  em = ups.children()[1];
	}
	else
	  continue;
	if(em.pid()<0) swap(ep,em);
      	LorentzTransform boostUps = LorentzTransform::mkFrameTransformFromBeta(p2.betaVec());
	FourMomentum pe = boost.transform(ep .momentum());
      	Vector3 axisE = boostUps.transform(pe).p3().unit();
      	axisPi.setZ(0.);
      	axisE.setZ(0.);
      	double chi = abs(atan2(axisE.cross(axisPi).dot(p2.p3().unit()), axisE.dot(axisPi)));
	if(chi>M_PI) chi=2.*M_PI-chi;
      	_h[1]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2010_I850492);

}

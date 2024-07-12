// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta -> pi+pi-pi0 and eta' -> pi0pi0pi0$
  class BESIII_2015_I1376484 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1376484);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==PID::ETA or
						Cuts::pid==PID::ETAPRIME);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      ETA.addStable(PID::K0S);
      ETA.addStable(PID::ETA);
      declare(ETA, "ETA");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h[ix  ],1,1,1+ix);
	book(_h[ix+2],2+ix,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1  = { {211,1}, {-211,1}, {111,1} };
      static const map<PdgId,unsigned int> & mode2   = { {111,3} };
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	// select right decay mode
	if ( ETA.decaying()[ix].pid()==PID::ETA && ETA.modeMatches(ix,3,mode1)) {
	  const Particle & pi0 = ETA.decayProducts()[ix].at( 111)[0];
	  const Particle & pip = ETA.decayProducts()[ix].at( 211)[0];
	  const Particle & pim = ETA.decayProducts()[ix].at(-211)[0];
	  double s1 = (pi0.momentum()+pim.momentum()).mass2();
	  double s2 = (pi0.momentum()+pip.momentum()).mass2();
	  double s3 = (pip.momentum()+pim.momentum()).mass2();
	  double mOut = pi0.mass()+pip.mass()+pim.mass();
	  double Q = ETA.decaying()[ix].mass()-mOut;
	  double X = sqrt(3.)/2./ETA.decaying()[ix].mass()/Q*(s1-s2);
	  double Y = 3.*(sqr(ETA.decaying()[ix].mass()-pi0.mass())-s3)/2./ETA.decaying()[ix].mass()/Q-1.;
	  _h[0]->fill(X);
	  _h[1]->fill(Y);
	}
	else if  ( ETA.modeMatches(ix,3,mode2)) {
	  const Particles & pi0 = ETA.decayProducts()[ix].at(111);
	  double s1 = (pi0[2].momentum()+pi0[1].momentum()).mass2();
	  double s2 = (pi0[2].momentum()+pi0[0].momentum()).mass2();
	  double s3 = (pi0[0].momentum()+pi0[1].momentum()).mass2();
	  double mOut = pi0[2].mass()+pi0[0].mass()+pi0[1].mass();
	  double Q = ETA.decaying()[ix].mass()-mOut;
	  double X = sqrt(3.)/2./ETA.decaying()[ix].mass()/Q*(s1-s2);
	  double Y = 3.*(sqr(ETA.decaying()[ix].mass()-pi0[2].mass())-s3)/2./ETA.decaying()[ix].mass()/Q-1.;
	  double Z = sqr(X)+sqr(Y);
	  if(ETA.decaying()[ix].pid()==PID::ETA) _h[2]->fill(Z);
	  else _h[3]->fill(Z);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1376484);

}

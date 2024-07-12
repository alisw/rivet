// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Jpsi -> gamma gamma phi
  class BESIII_2018_I1646687 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1646687);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::PHI);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      book(_h_mass,1,1,1);
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_angle[ix],2,1,1+ix);
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
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode = { { 22,2},{ 333,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if( PSI.decaying().size()!=1) vetoEvent;
      if(!PSI.modeMatches(0,3,mode)) vetoEvent;
      // particles
      const Particle  & phi = PSI.decayProducts()[0].at(333)[0];
      const Particles & gam = PSI.decayProducts()[0].at( 22);
      for(unsigned int ix=0;ix<2;++ix) {
	double mPhiGamma = (gam[ix].momentum()+phi.momentum()).mass();
	_h_mass->fill(mPhiGamma);
	double cTheta = axis.dot(gam[ix==0 ? 1 : 0].p3().unit());
	if(mPhiGamma>1.4 && mPhiGamma<1.6)
	  _h_angle[0]->fill(cTheta);
	else if(mPhiGamma>1.75 && mPhiGamma<1.9)
	  _h_angle[1]->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass,1.,false);
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_angle[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass, _h_angle[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1646687);

}

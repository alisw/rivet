// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c -> phi K K pi
  class BESIII_2015_I1352828 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1352828);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==10441 or
						Cuts::pid==20443 or
						Cuts::pid==445);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PHI);
      chi.addStable( PID::K0S);
      chi.addStable( PID::PI0);
      declare(chi, "chi");
      // histos
      for(unsigned int iy=0;iy<2;++iy)
	for(unsigned int ix=0;ix<3;++ix)
	  book(_h[iy][ix],1,1+iy,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 333,1}, { 310,1}, { 321,1}, {-211,1} };
      static const map<PdgId,unsigned int> & mode1CC = { { 333,1}, { 310,1}, {-321,1}, { 211,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 333,1}, { 321,1}, {-321,1}, { 111,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	int imode = 0;
	if(chi.modeMatches(ix,4,mode1) ||
	   chi.modeMatches(ix,4,mode1CC))    imode=0;
	else if(chi.modeMatches(ix,4,mode2)) imode=1;
	else continue;
	unsigned int ichi = 0;
	if     (chi.decaying()[ix].pid()==20443) ichi=1;
	else if(chi.decaying()[ix].pid()==  445) ichi=2;
	const Particle & phi = chi.decayProducts()[ix].at(333)[0];
	double mass = (chi.decaying()[ix].momentum()-phi.momentum()).mass();
	_h[imode][ichi]->fill(mass);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=0;iy<2;++iy)
	for(unsigned int ix=0;ix<3;++ix)
	  normalize(_h[iy][ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1352828);

}

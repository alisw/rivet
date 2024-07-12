// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> pi,rho l+ nu
  class BABAR_2011_I855306 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2011_I855306);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==PID::BPLUS ||
				Cuts::abspid==PID::B0   ), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_c[ix],"TMP/c_"+toString(ix+1));
	for(unsigned int iy=0;iy<3;++iy) {
	  book(_h[ix][iy],ix+1,1,iy+1);
	}
      }
    }
    
    // Calculate the Q2 using mother and daughter charged lepton
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - filter_select(B.children(), Cuts::pid==mesonID)[0];
      return q*q;
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get B+ Mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	// B0 modes
	if(p.abspid()==511) {
	  _c[0]->fill();
	  if (isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E}) ||
	      isSemileptonicDecay(p, {PID::PIMINUS, PID::ANTIMUON, PID::NU_MU}) ) {
	    double qq = q2(p,PID::PIMINUS);
            _h[0][0]->fill(qq);
            _h[0][2]->fill(qq);
	  }
	  else if(isSemileptonicDecay(p, {PID::PIPLUS , PID::ELECTRON, PID::NU_EBAR}) ||
		  isSemileptonicDecay(p, {PID::PIPLUS , PID::MUON    , PID::NU_MUBAR})) {
	    double qq = q2(p,PID::PIPLUS);
            _h[0][0]->fill(qq);
            _h[0][2]->fill(qq);
	  }
	  else if (isSemileptonicDecay(p, {PID::RHOMINUS, PID::POSITRON, PID::NU_E}) ||
		   isSemileptonicDecay(p, {PID::RHOMINUS, PID::ANTIMUON, PID::NU_MU})) {
	    double qq = q2(p,PID::RHOMINUS);
            _h[1][0]->fill(qq);
            _h[1][2]->fill(qq);
	  }
	  else if( isSemileptonicDecay(p, {PID::RHOPLUS, PID::ELECTRON, PID::NU_EBAR}) ||
		   isSemileptonicDecay(p, {PID::RHOPLUS, PID::MUON    , PID::NU_MUBAR})) {
	    double qq = q2(p,PID::RHOPLUS);
            _h[1][0]->fill(qq);
            _h[1][2]->fill(qq);
	  }
	}
	// B+ modes
	else {
	  _c[1]->fill();
	  if (isSemileptonicDecay(p, {PID::PI0, PID::POSITRON, PID::NU_E}) ||
	      isSemileptonicDecay(p, {PID::PI0, PID::ANTIMUON, PID::NU_MU}) ||
	      isSemileptonicDecay(p, {PID::PI0, PID::ELECTRON, PID::NU_EBAR}) ||
	      isSemileptonicDecay(p, {PID::PI0, PID::MUON    , PID::NU_MUBAR})) {
            _h[0][1]->fill(q2(p,PID::PI0));
	  }
	  else if (isSemileptonicDecay(p, {PID::RHO0, PID::POSITRON, PID::NU_E}) ||
		   isSemileptonicDecay(p, {PID::RHO0, PID::ANTIMUON, PID::NU_MU}) ||
		   isSemileptonicDecay(p, {PID::RHO0, PID::ELECTRON, PID::NU_EBAR}) ||
		   isSemileptonicDecay(p, {PID::RHO0, PID::MUON    , PID::NU_MUBAR})) {
            _h[1][1]->fill(q2(p,PID::RHO0));
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Br in units 10^-4 and average over e/mu modes
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  if(iy==1) scale(_h[ix][iy], 0.5*1e4/ *_c[1]);
	  else      scale(_h[ix][iy], 0.5*1e4/ *_c[0]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2011_I855306);

}

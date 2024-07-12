// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"


namespace Rivet {


  /// @brief J/psi psi(2S) -> p pbar eta'
  class BESIII_2019_I1709205 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1709205);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443 or Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<6;++iy) {
	  book(_h[ix][iy],1+ix,1,1+iy);
	}
      }
      book(_dalitz[0], "dalitz_1",50,3.,8.,50,3.,8.);
      book(_dalitz[1], "dalitz_2",50,3.5,5.,50,3.5,5.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 2212,1}, {-2212,1}, { 331,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particles & etap = psi.decayProducts()[ix].at( 331);
	const Particles & pp   = psi.decayProducts()[ix].at( 2212);
	const Particles & pbar = psi.decayProducts()[ix].at(-2212);
	double mminus = (pbar[0].momentum()+etap[0].momentum()).mass2();
	double mplus  = (pp  [0].momentum()+etap[0].momentum()).mass2();
	double mneut  = (pp  [0].momentum()+pbar[0].momentum()).mass2();
	unsigned int iloc = psi.decaying()[ix].pid()==443 ? 1 : 0;
	for(unsigned int ix=0;ix<4;ix+=3) {
	  _h[iloc][ix+2]->fill(sqrt(mneut ));
	  _h[iloc][ix+0]->fill(sqrt(mplus ));
	  _h[iloc][ix+1]->fill(sqrt(mminus));
	}
	_dalitz[iloc]->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_dalitz[ix]);
	for(unsigned int iy=0;iy<6;++iy) {
	  normalize(_h[ix][iy]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][6];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1709205);

}

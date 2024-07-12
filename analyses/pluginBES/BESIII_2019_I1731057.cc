// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> K+K-pi0
  class BESIII_2019_I1731057 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1731057);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::PI0);
      declare(PSI,"PSI");
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_KpKm [ix],1+ix,1,1);
	book(_h_Kppi0[ix],1+ix,1,2);
      }
      book(_dalitz, "dalitz",50,0.,7.,50,0.0,7.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-321,1}, {111,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for(unsigned int ix=0;ix<PSI.decaying().size();++ix) {
	if (!PSI.modeMatches(ix,3,mode)) continue;
	const Particle & Kp  = PSI.decayProducts()[ix].at( 321)[0];
	const Particle & Km  = PSI.decayProducts()[ix].at(-321)[0];
	const Particle & pi0 = PSI.decayProducts()[ix].at( 111)[0];
	double mminus = (Km.momentum()+pi0.momentum()).mass2();
	double mplus  = (Kp.momentum()+pi0.momentum()).mass2();
	double mneut  = (Kp.momentum()+Km.momentum()).mass2();
	_dalitz->fill(mplus,mminus);
	mplus=sqrt(mplus);
	mminus=sqrt(mminus);
	mneut=sqrt(mneut);
	_h_KpKm [0]->fill(mneut );
	_h_Kppi0[0]->fill(mplus );
	_h_Kppi0[0]->fill(mminus);
	if(mplus>1.05 && mminus>1.05) {
	  _h_KpKm [1]->fill(mneut );
	  _h_Kppi0[1]->fill(mplus );
	  _h_Kppi0[1]->fill(mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_KpKm [ix],1.,false);
	normalize(_h_Kppi0[ix],1.,false);
      }
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KpKm[2],_h_Kppi0[2];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1731057);

}

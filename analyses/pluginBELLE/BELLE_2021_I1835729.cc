// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Xi_c0 -> Xi0 K+K-
  class BELLE_2021_I1835729 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1835729);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4132);
      declare(ufs, "UFS");
      DecayedParticles XIC0(ufs);
      XIC0.addStable(PID::PI0);
      XIC0.addStable(PID::K0S);
      XIC0.addStable(PID::ETA);
      XIC0.addStable(PID::ETAPRIME);
      XIC0.addStable(PID::XI0);
      declare(XIC0, "XIC0");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
      book(_dalitz, "dalitz",50,0.9,1.4,50,3.2,4.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {-321,1},{ 321,1}, { PID::XI0,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1},{-321,1}, {-PID::XI0,1}};
      DecayedParticles XIC0 = apply<DecayedParticles>(event, "XIC0");
      // loop over particles
      for(unsigned int ix=0;ix<XIC0.decaying().size();++ix) {
	int sign = 1;
	if (XIC0.decaying()[ix].pid()>0 && XIC0.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (XIC0.decaying()[ix].pid()<0 && XIC0.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & xi0 = XIC0.decayProducts()[ix].at( sign*PID::XI0);
	const Particles & Kp  = XIC0.decayProducts()[ix].at( sign*PID::KPLUS);
	const Particles & Km  = XIC0.decayProducts()[ix].at( sign*PID::KMINUS);
	double mXiKp = (xi0[0].momentum()+Kp[0].momentum()).mass2();
	double mXiKm = (xi0[0].momentum()+Km[0].momentum()).mass2();
	double mKK   = (Kp [0].momentum()+Km[0].momentum()).mass2();
	_dalitz ->fill(mKK,mXiKm);
	_h[0]->fill(sqrt(mKK  ));
	_h[1]->fill(sqrt(mXiKm));
	_h[2]->fill(sqrt(mXiKp));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix]);
      normalize(_dalitz );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1835729);

}

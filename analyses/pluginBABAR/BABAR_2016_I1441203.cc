// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> pi+pi-pi0
  class BABAR_2016_I1441203 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2016_I1441203);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // histograms
      book(_h_pi[0],1,1,1);
      book(_h_pi[1],1,1,2);
      book(_h_pi[2],1,1,3);
      book(_dalitz, "dalitz",50,0.,3.2,50,0.0,3.2);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{-211,1}, {111,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if( ! D0.modeMatches(ix,3,mode)  ) continue;
	int sign = D0.decaying()[ix].pid()/421;
	const Particles & pi0 = D0.decayProducts()[ix].at(111);
	const Particles & pip = D0.decayProducts()[ix].at( sign*211);
	const Particles & pim = D0.decayProducts()[ix].at(-sign*211);
	double mneut = (pim[0].momentum()+pip[0].momentum()).mass2();
	// if(mneut>475*MeV && mneut<505*MeV) continue;
	double mplus  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	double mminus = (pim[0].momentum()+pi0[0].momentum()).mass2();
	_h_pi[0]->fill(mplus);
	_h_pi[1]->fill(mminus);
	_h_pi[2]->fill(mneut);
	_dalitz->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h_pi[ix]);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi[3];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2016_I1441203);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> eta eta pi0
  class BESIII_2018_I1662665 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1662665);


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
      book(_h_pieta,1,1,1);
      book(_dalitz,"dalitz",50.,0.4,1.6,50,0.4,1.8);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { { 221,2}, {111,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if ( ! D0.modeMatches(ix,3,mode)) continue;
	const Particles & eta = D0.decayProducts()[ix].at(221);
	const Particles & pi0 = D0.decayProducts()[ix].at(111);
	double m1 = (eta[0].momentum()+pi0[0].momentum()).mass2();
	double m2 = (eta[1].momentum()+pi0[0].momentum()).mass2();
	if(m1>m2) swap(m1,m2);
	_dalitz->fill(m1,m2);
	_h_pieta->fill(sqrt(m1));
	_h_pieta->fill(sqrt(m2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pieta);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pieta;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1662665);

}

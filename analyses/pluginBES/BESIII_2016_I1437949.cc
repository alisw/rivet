// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma eta' pi+pi-
  class BESIII_2016_I1437949 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2016_I1437949);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::ETA);
      PSI.addStable(PID::ETAPRIME);
      declare(PSI, "PSI");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
      book(_h[2],2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode = { { 22,1}, { 331,1}, {211,1}, {-211,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      for(unsigned int ix=0;ix<PSI.decaying().size();++ix) {
	if(!PSI.modeMatches(ix,4,mode)) continue;
	const Particle  & gam  = PSI.decayProducts()[ix].at( 22)[0];
	double mass =(PSI.decaying()[ix].momentum()-gam.momentum()).mass() ;
	for(unsigned int ix=0;ix<3;++ix) _h[ix]->fill(mass);
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


  RIVET_DECLARE_PLUGIN(BESIII_2016_I1437949);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta' -> pi+pi- e+e- decay
  class BESIII_2020_I1830421 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1830421);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 331);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      declare(ETA,"ETA");
      // Book histograms
      //book(_h_mee  , 1, 1, 1);
      book(_h_mpipi, 1, 1, 2);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 211,1}, {-211,1}, { 11,1}, {-11,1}};
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	  if (!ETA.modeMatches(ix,4,mode)) continue; 
	  const Particle & pip = ETA.decayProducts()[ix].at( 211)[0];
	  const Particle & pim = ETA.decayProducts()[ix].at(-211)[0];
	  //const Particle & ep  = ETA.decayProducts()[ix].at(  11)[0];
	  //const Particle & em  = ETA.decayProducts()[ix].at(- 11)[0];
	  //_h_mee  ->fill((em .momentum()+ep .momentum()).mass());
	  _h_mpipi->fill((pim.momentum()+pip.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //normalize(_h_mee);
      normalize(_h_mpipi);
    }

    /// @}


    /// @name Histograms
    /// @{
    //Histo1DPtr _h_mee;
    Histo1DPtr _h_mpipi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1830421);

}

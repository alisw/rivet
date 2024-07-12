// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Ds+ -> K+ pi+ pi-
  class BESIII_2022_I2084294 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2084294);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      DS.addStable(PID::ETA);
      DS.addStable(PID::ETAPRIME);
      declare(DS, "DS");
      // histograms
      book(_h_Kpip,1,1,1);
      book(_h_Kpim,1,1,2);
      book(_h_pipi,1,1,3);
      book(_dalitz,"dalitz",50,0.,2.3,50,0.3,3.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1}, { 211,1},{-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1}, {-211,1},{ 211,1}};
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for(unsigned int ix=0;ix<DS.decaying().size();++ix) {
	int sign = 1;
	if (DS.decaying()[ix].pid()>0 && DS.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (DS.decaying()[ix].pid()<0 && DS.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & Kp = DS.decayProducts()[ix].at( sign*321)[0];
	const Particle & pip= DS.decayProducts()[ix].at( sign*211)[0];
	const Particle & pim= DS.decayProducts()[ix].at(-sign*211)[0];
	double mp    = (Kp .momentum()+pip.momentum()).mass2();
	double mm    = (Kp .momentum()+pim.momentum()).mass2();
	double mpipi = (pip.momentum()+pim.momentum()).mass2();
	_dalitz->fill(mpipi,mm);
	if(mpipi<sqr(.4676) || mpipi>sqr(0.5276)) {
	  _h_pipi->fill(sqrt(mpipi));
	  _h_Kpip->fill(sqrt(mp));
	  _h_Kpim->fill(sqrt(mm));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpip);
      normalize(_h_Kpim);
      normalize(_h_pipi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpip,_h_Kpim,_h_pipi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2084294);

}

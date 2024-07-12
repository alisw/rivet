// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D0 -> K- pi+ pi+ pi-
  class BESIII_2017_I1511280 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1511280);


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
      for(unsigned int ix=0;ix<8;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { {-321,1},{ 211,2}, {-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1},{-211,2}, { 211,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = 1;
	if (D0.decaying()[ix].pid()>0 && D0.modeMatches(ix,4,mode)) {
	  sign=1;
	}
	else if  (D0.decaying()[ix].pid()<0 && D0.modeMatches(ix,4,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particles & Km = D0.decayProducts()[ix].at(-sign*321);
	const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
	// first pip gives higher pi+pi- mass
	double mpipi[2] = {(pip[0].momentum()+pim[0].momentum()).mass(),
			   (pip[1].momentum()+pim[0].momentum()).mass()};
	unsigned int i1= mpipi[0]>mpipi[1] ? 0 : 1;
	unsigned int i2= i1==0 ? 1 : 0;
	_h[0]->fill((Km [0].momentum()+pip[i1].momentum()).mass());
	_h[1]->fill((Km [0].momentum()+pip[i2].momentum()).mass());
	_h[2]->fill(sqr(mpipi[i1]));
	_h[3]->fill(sqr(mpipi[i2]));
	_h[4]->fill((Km [0].momentum()+pim[0].momentum()+pip[i1].momentum()).mass());
	_h[5]->fill((Km [0].momentum()+pim[0].momentum()+pip[i2].momentum()).mass());
	_h[6]->fill((pim[0].momentum()+pip[i1].momentum()+pip[i2].momentum()).mass());
	_h[7]->fill((Km [0].momentum()+pip[i1].momentum()+pip[i2].momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<8;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[8];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1511280);

}

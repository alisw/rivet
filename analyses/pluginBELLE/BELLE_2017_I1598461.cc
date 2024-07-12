// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> K+ K- pi+
  class BELLE_2017_I1598461 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2017_I1598461);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      declare(BP, "BP");
      // book histograms
      for(unsigned int ix=0;ix<3;++ix) {
	if(ix==0) book(_h[ix],1,1,1);
	else      book(_h[ix],"TMP/h_"+toString(ix+1),refData(1,1,1));
	book(_c[ix],"TMP/c_"+toString(ix+1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-321,1}, { 211,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 321,1},{ 321,1}, {-211,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	_c[0]->fill();
	if(BP.decaying()[ix].pid()>0) _c[1]->fill();
	else                          _c[2]->fill();
      	int sign = 1;
      	if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
      	const Particle & Kp  = BP.decayProducts()[ix].at( 321)[0];
      	const Particle & Km  = BP.decayProducts()[ix].at(-321)[0];
	double mKK = (Kp.momentum()+Km.momentum()).mass();
	_h[0]->fill(mKK);
	if(sign==1) _h[1]->fill(mKK);
	else        _h[2]->fill(mKK);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	scale(_h[ix],1e7/ *_c[ix]);
      Scatter2DPtr as;
      book(as,1,1,2);
      asymm(_h[2],_h[1],as);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2017_I1598461);

}

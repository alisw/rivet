// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief  B-> K pi psi(2S)
  class BELLE_2009_I820737 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I820737);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 ||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(310);
      BB.addStable(100443);
      declare(BB, "BB");
      // histograms
      vector<double> bins1 = {0,sqr(0.796),sqr(0.996),sqr(1.332),sqr(1.532),3};
      vector<double> bins2 = {0,19.0,20.5,23};
      for(unsigned int ix=0;ix<5;++ix) {
	Histo1DPtr tmp;
	book(tmp,1,1,1+ix);
	_h_piPsi.add(bins1[ix],bins1[ix+1],tmp);
	if(ix>2) continue;
	book(tmp,2,1,1+ix);
	_h_Kpi.add(bins2[ix],bins2[ix+1],tmp);
      }
      book(_h_veto,3,1,1);
      book(_c,"TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,1},{-211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-321,1},{ 211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 310,1},{-211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 310,1},{ 211,1}, { 100443,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int sign = 1,iK(0);
	if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
	  sign=1; iK = 321;
	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
       	  sign=-1; iK=-321;
       	}
	else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
	  sign=1; iK = 310;
	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
       	  sign=-1; iK= 310;
       	}
      	else
      	  continue;
	_c->fill();
       	const Particle & Kp  = BB.decayProducts()[ix].at( iK)[0];
       	const Particle & pim = BB.decayProducts()[ix].at(-211*sign)[0];
       	const Particle & psi = BB.decayProducts()[ix].at( 100443  )[0];
       	double m2Kpi  = (Kp .momentum()+pim.momentum()).mass2();
      	double m2Psipi= (psi.momentum()+pim.momentum()).mass2();
	_h_piPsi.fill(m2Kpi  ,m2Psipi);
	_h_Kpi  .fill(m2Psipi,m2Kpi  );
	if(m2Kpi<sqr(0.796) || (m2Kpi>sqr(0.996)&&m2Kpi<sqr(1.332)) ||
	   m2Kpi>sqr(1.532)) _h_veto->fill(m2Psipi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(Histo1DPtr h : _h_piPsi.histos()) {
	scale(h, 1./ *_c);
      }
      for(Histo1DPtr h : _h_Kpi.histos()) {
	scale(h, 1./ *_c);
      }
      normalize(_h_veto,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_piPsi,_h_Kpi;
    CounterPtr _c;
    Histo1DPtr _h_veto;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I820737);

}

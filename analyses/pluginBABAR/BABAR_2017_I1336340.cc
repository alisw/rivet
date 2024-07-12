// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> KS0 pi+ pi0
  class BABAR_2017_I1336340 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2017_I1336340);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(310);
      BP.addStable(111);
      declare(BP, "BP");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_sum[ix],1,1,1+ix);
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h_charge[iy][ix],2,1+iy,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, { 111,1}, { 310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1}, { 111,1}, { 310,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign = 1;
      	if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode))        sign= 1;
      	else if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) sign=-1;
      	else continue;
      	const Particle & pip = BP.decayProducts()[ix].at( sign*211)[0];
  	const Particle & pi0 = BP.decayProducts()[ix].at(      111)[0];
  	const Particle & K0  = BP.decayProducts()[ix].at(      310)[0];
	double mKpip = (K0 .momentum()+pip.momentum()).mass();
	double mKpi0 = (K0 .momentum()+pi0.momentum()).mass();
	// veto D decaysa
	if(mKpi0>1.804 && mKpi0<1.924) continue;
	double mpipi = (pi0.momentum()+pip.momentum()).mass();
	_h_sum[0]->fill(mKpip);
	_h_sum[1]->fill(mKpi0);
	_h_sum[2]->fill(mpipi);
	_h_charge[(1-sign)/2][0]->fill(mKpip);
	_h_charge[(1-sign)/2][1]->fill(mKpi0);
	_h_charge[(1-sign)/2][2]->fill(mpipi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_h_sum[ix],1.,false);
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h_charge[iy][ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sum[3],_h_charge[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2017_I1336340);

}

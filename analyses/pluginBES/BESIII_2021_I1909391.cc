// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D_s+ -> pi+pi+pi-
  class BESIII_2021_I1909391 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1909391);


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
      declare(DS,"DS");
      // histos
      book(_h_pippim[0],1,1,1);
      book(_h_pippip   ,1,1,2);
      book(_h_pippim[1],1,1,3);
      book(_h_pippim[2],1,1,4);
      book(_dalitz, "dalitz",50,0.,3.5,50,0.0,3.5);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // parameters from the efficiency function, table 1 in paper
      static const double E1  =  0.064;
      static const double E2  = -0.066;
      static const double E3  = -0.006;
      static const double E11 = -0.158;
      static const double E12 =  0.090;
      static const double Eth[3] = {1.516,1.516,1.563};
      static const map<PdgId,unsigned int> & mode   = { { 211,2}, {-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,2}, { 211,1}};
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
	const Particles & pip = DS.decayProducts()[ix].at( sign*211);
	const Particle  & pim = DS.decayProducts()[ix].at(-sign*211)[0];
	// kinematic variables
	double x[3] = {(pim.momentum()+pip[0].momentum()).mass2(),
		       (pim.momentum()+pip[1].momentum()).mass2(),
		       (pip[0].momentum()+pip[1].momentum()).mass2()};
	if(x[0]>x[1]) swap(x[0],x[1]);
	_dalitz->fill(x[0],x[1]);
	_dalitz->fill(x[1],x[0]);
	// calculate the efficiency
	double xh = x[0]-1.,yh=x[1]-1.;
	double eff = (1.+E1*(xh+yh)+E2*(sqr(xh)+sqr(yh))+E3*(pow(xh,3)+pow(yh,3))
		      +E11*xh*yh+E12*xh*yh*(xh+yh));
	double xmax = sqr(DS.decaying()[ix].mass()-pip[0].mass());
	double T=1.;
	for(unsigned int ix=0;ix<3;++ix) {
	  double arg = Eth[ix]*abs(x[ix]-xmax);
	  if(arg<0.5*M_PI) T *=sin(arg);
	}
	eff *=T;
	// fill plots
	_h_pippim[0]->fill(x[0],eff);
	_h_pippim[0]->fill(x[1],eff);
	_h_pippim[1]->fill(x[0],eff);
	_h_pippim[2]->fill(x[1],eff);
	_h_pippip   ->fill(x[2],eff);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h_pippim[ix]);
      normalize(_h_pippip);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim[3],_h_pippip;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1909391);

}

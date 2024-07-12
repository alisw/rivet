// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> K- pi+ pi+
  class CLEOC_2008_I780363 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOC_2008_I780363);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles DP(ufs);
      DP.addStable(PID::PI0);
      DP.addStable(PID::K0S);
      DP.addStable(PID::ETA);
      DP.addStable(PID::ETAPRIME);
      declare(DP, "DP");
      // histos
      book(_h_Kpiall,1,1,1);
      book(_h_pipi  ,1,2,1);

    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // parameters from the efficiency function, table 1 in paper
      static const double E1   = -0.0153;
      static const double E2   = -0.030;
      static const double E3   =  0.162;
      static const double Exy  = -0.053;
      static const double Exyn =  0.673;
      // static const double Eth[3] = {4.25,4.25,2.907};


      static const map<PdgId,unsigned int> & mode   = { { 211,2},{-321,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,2},{ 321,1}};
      DecayedParticles DP = apply<DecayedParticles>(event, "DP");
      // loop over particles
      for(unsigned int ix=0;ix<DP.decaying().size();++ix) {
	int sign = 1;
	if (DP.decaying()[ix].pid()>0 && DP.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (DP.decaying()[ix].pid()<0 && DP.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle  & Km = DP.decayProducts()[ix].at(-sign*321)[0];
	const Particles & pip= DP.decayProducts()[ix].at( sign*211);
	// kinematic variables
	double x[3]  = {(Km    .momentum() +pip[0].momentum()).mass2(),
			(Km    .momentum() +pip[1].momentum()).mass2(),
			(pip[0].momentum()+pip[1].momentum()).mass2()};
	if(x[1]<x[0]) swap(x[0],x[1]);
	// calculate the efficiency, eqns 6,7 from paper
	// double xmax[3] = {sqr(meson.mass()-pip[0].mass()),
	// 		    sqr(meson.mass()-pip[0].mass()),
	// 		    sqr(meson.mass()-Km  .mass())};
	double xh = x[0]-1.5,yh=x[1]-1.5;
	double eff = (1.+E1*(xh+yh)+E2*(sqr(xh)+sqr(yh))+E3*(pow(xh,3)+pow(yh,3))
		      +Exy*xh*yh+Exyn*xh*yh*(xh+yh));
	// double T=1.;
	// for(unsigned int ix=2;ix<3;++ix) {
	//   double arg = Eth[ix]*abs(x[ix]-xmax[ix]);
	//   if(arg<0.5*M_PI) T *=sin(arg);
	// }
	// eff *=T;
	// fill plots
	_h_Kpiall->fill( x[0],eff);
	_h_Kpiall->fill( x[1],eff);
	_h_pipi  ->fill( x[2],eff);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpiall);
      normalize(_h_pipi  );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Kpiall, _h_pipi;
    //@}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2008_I780363);

}

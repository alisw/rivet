// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s -> pi+ pi+ pi- +X
  class BESIII_2022_I2618227 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2618227);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==431),"UFS");
      // histos
      book(_h_br,1,1,1);
      book(_h_mass,2,1,1);
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_mom[ix],3,1,1+ix);
      book(_c,"TMP/nDs");
    }
    
    void findChildren(const Particle & p, Particles & piPlus, Particles & piMinus) {
      for( const Particle &child : p.children()) {
	if(child.pid()==PID::PIPLUS) {
	  piPlus.push_back(child);
	}
	else if(child.pid()==PID::PIMINUS) {
	  piMinus.push_back(child);
	}
	else if(child.pid()==PID::K0S || child.pid()==PID::K0L)
	  continue;
	else if(!child.children().empty())
	  findChildren(child,piPlus,piMinus);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (const Particle& p :  apply<FinalState>(event, "UFS").particles()) {
	_c->fill();
	Particles piPlus,piMinus;
	findChildren(p,piPlus,piMinus);
	if(p.pid()<0) swap(piPlus,piMinus);
	if(piPlus.size()>=2 && !piMinus.empty()) {
	  _h_br->fill(.5);
	  sortBy(piMinus,cmpMomByP);
	  sortBy(piPlus ,cmpMomByP);
	  _h_mom[0]->fill(piMinus[0].p3().mod());
	  for(unsigned int ix=0;ix<2;++ix)
	    _h_mom[1+ix]->fill(piPlus[ix].p3().mod());
	  double mass = (piMinus[0].momentum()+piPlus[0].momentum()+piPlus[1].momentum()).mass();
	  _h_mass->fill(mass);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_br,   100./ *_c);
      scale(_h_mass, 100./ *_c);
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h_mom[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_br,_h_mass,_h_mom[3];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2618227);

}

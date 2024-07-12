// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0,+ -> pi+ pi+ pi- +X
  class BESIII_2023_I2621481 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2621481);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==411 or
				Cuts::abspid==421),"UFS");
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_br  [ix],1,1,1+ix);
	book(_h_mass[ix],2,1,1+ix);
	for(unsigned int iy=0;iy<3;++iy)
	  book(_h_mom[ix][iy],3,1+ix,1+iy);
	book(_c[ix],"TMP/nD_"+toString(ix+1));
      }
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
	unsigned int iloc = p.abspid()==421 ? 0 : 1;
	_c[iloc]->fill();
	Particles piPlus,piMinus;
	findChildren(p,piPlus,piMinus);
	if(p.pid()<0) swap(piPlus,piMinus);
	if(piPlus.size()>=2 && !piMinus.empty()) {
	  _h_br[iloc]->fill(.5);
	  sortBy(piMinus,cmpMomByP);
	  sortBy(piPlus ,cmpMomByP);
	  _h_mom[iloc][2]->fill(piMinus[0].p3().mod());
	  for(unsigned int ix=0;ix<2;++ix)
	    _h_mom[iloc][ix]->fill(piPlus[ix].p3().mod());
	  double mass = (piMinus[0].momentum()+piPlus[0].momentum()+piPlus[1].momentum()).mass();
	  _h_mass[iloc]->fill(mass);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_br  [ix], 100./ *_c[ix]);
	scale(_h_mass[ix], 100./ *_c[ix]);
	for(unsigned int iy=0;iy<3;++iy)
	  normalize(_h_mom[ix][iy],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_br[2],_h_mass[2],_h_mom[2][3];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2621481);

}

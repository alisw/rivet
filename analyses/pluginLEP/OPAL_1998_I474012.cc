// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief b-baryon polarization
  class OPAL_1998_I474012 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_1998_I474012);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_El   , "El",   45,0.,45.0);
      book(_h_Ev   , "Ev",   45,0.,45.0);
      book(_h_ratio, "ratio",20,0.,10.0);

    }

    void findDecayProducts(Particle p, Particles & lep, Particles & nu) {
      for(const Particle & child : p.children()) {
	if(PID::isHadron(child.pid())) continue;
	if(child.abspid()==11 or child.abspid()==13)
	  lep.push_back(child);
	else if(child.abspid()==12 or child.abspid()==14)
	  nu.push_back(child);
	else if(child.abspid()!=15)
	  findDecayProducts(child,lep,nu);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over weakly decaying b-baryons
      for (const Particle& p : ufs.particles(Cuts::abspid==5122)) {
	Particles lep,nu;
	findDecayProducts(p,lep,nu);
	if(lep.size()!=1 || nu.size()!=1) continue;
	_h_El   ->fill(lep[0].momentum().t());
	_h_Ev   ->fill(nu [0].momentum().t());
	_h_ratio->fill(nu [0].momentum().t()/lep[0].momentum().t());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_El   );
      normalize(_h_Ev   );
      normalize(_h_ratio);
      if(_h_El->effNumEntries()!=0. and _h_Ev->effNumEntries()!=0.) {
	double Ev  = _h_Ev->xMean();
	double El  = _h_El->xMean();
	double dEv = _h_Ev->xStdErr();
	double dEl = _h_El->xStdErr();
	double ratio = Ev/El;
	double dr    = (El*dEv-Ev*dEl)/sqr(El);
	double rho = 0.091;
	double  P = 7. - 20./(2. + ratio) + 10.*rho*ratio*(-4. + 3.*ratio)/sqr(2.+ratio);
	double dP = (20.*(2. + ratio + rho*(-4. + 8.*ratio)))/pow(2.+ratio,3)*dr;
       	Scatter2DPtr h_pol;
	book(h_pol, 1,1,1);
       	h_pol->addPoint(91.2, P, make_pair(0.5,0.5),make_pair(dP,dP) );
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_El, _h_Ev, _h_ratio;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(OPAL_1998_I474012);


}

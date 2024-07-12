// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief phi and D_s spectra at 29 GeV
  class HRS_1985_I213242 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1985_I213242);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_phi,1,1,1);
      book(_h_Ds ,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for(const Particle & p : ufs.particles(Cuts::abspid==431 ||
					     Cuts::pid==333)) {
      	double xE = 2.*p.E()/sqrtS();
      	Vector3 mom3 = p.p3();
        const double energy = p.E();
      	double modp = mom3.mod();
      	double beta = modp/energy;
	if(p.pid()==333) 
	  _h_phi->fill(xE,1./beta);
	else
	  _h_Ds ->fill(xE,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // PDG 2020 D_s phi pi br
      double br=0.045;
      scale( _h_phi, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale( _h_Ds , br*sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_phi,_h_Ds;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1985_I213242);

}

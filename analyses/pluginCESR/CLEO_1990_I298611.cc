// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_c spectra
  class CLEO_1990_I298611 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1990_I298611);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<6;++ix) {
	_h_x.push_back(Histo1DPtr());
	book(_h_x[ix],ix+1,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.28646));
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==4122)) {
	double xp = p.momentum().p3().mod()/Pmax;
        for(unsigned int ix=0;ix<6;++ix)
	  _h_x[ix]->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // branching ratios for the different modes used in the measurements
      double br[6]={0.0623,2.*0.0158,2.*0.0159,
		    0.0129,0.0361   ,0.0062};
      for(unsigned int ix=0;ix<6;++ix)
	scale(_h_x[ix],crossSection()*br[ix]/picobarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    vector<Histo1DPtr> _h_x;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1990_I298611);

}

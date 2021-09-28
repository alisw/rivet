// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief baryons at 34.8 and 42.1 GeV
  class TASSO_1989_I266893 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1989_I266893);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Sphericity(cfs), "Sphericity");
      // Book histograms
      _ih=-1;
      if(fuzzyEquals(sqrtS()/GeV, 34.8, 1e-3)) {
	_ih=0;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 42.1, 1e-3)) {
	_ih=1;
      }
      else
	MSG_ERROR("Beam energy " << sqrtS() << " not supported!");

      book(_h_lam_p    ,6*_ih+3,1,1);
      book(_h_lam_pL   ,6*_ih+4,1,1);
      book(_h_lam_pTIn ,6*_ih+5,1,1);
      book(_h_lam_pTOut,6*_ih+6,1,1);
      book(_h_lam_rap  ,6*_ih+7,1,1);
      book(_h_lam_x    ,6*_ih+8,1,1);
      book(_p_lam_S_1  ,15+_ih,1,1);
      book(_p_lam_S_2  ,15+_ih,1,2);
      if(_ih==0) {
      	book(_h_xi_p    ,18,1,1);
      	book(_h_xi_pL   ,19,1,1);
      	book(_h_xi_pTIn ,20,1,1);
      	book(_h_xi_pTOut,21,1,1);
      	book(_h_xi_rap  ,22,1,1);
      	book(_h_xi_x    ,23,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t numParticles = cfs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
      				   beams.second.p3().mod() ) / 2.0;
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      unsigned int nLam(0);
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for(const Particle & p : ufs.particles(Cuts::abspid==3122 or Cuts::abspid==3312)) {
      	int id = abs(p.pid());
      	double xE = p.E()/meanBeamMom;
      	Vector3 mom3 = p.p3();
        const double energy = p.E();
      	double modp = mom3.mod();
      	double beta = modp/energy;
        const double momS = dot(sphericity.sphericityAxis(), mom3);
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
        const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));
      	if(id==3122) {
      	  _h_lam_x->fill(xE,1./beta);
      	  _h_lam_p->fill(modp/GeV);
      	  _h_lam_pL   ->fill(abs(momS)/GeV  );
      	  _h_lam_pTIn ->fill(abs(pTinS)/GeV );
      	  _h_lam_pTOut->fill(abs(pToutS)/GeV);
      	  _h_lam_rap  ->fill(abs(rapidityS) );
	  ++nLam;
      	}
      	else if(_h_xi_x) {
      	  _h_xi_x->fill(xE,1./beta);
      	  _h_xi_p->fill(modp/GeV);
      	  _h_xi_pL   ->fill(abs(momS)/GeV  );
      	  _h_xi_pTIn ->fill(abs(pTinS)/GeV );
      	  _h_xi_pTOut->fill(abs(pToutS)/GeV);
      	  _h_xi_rap  ->fill(abs(rapidityS) );
      	}
      }
      double sphere = sphericity.sphericity();
      _p_lam_S_1->fill(sphere,nLam);
      _p_lam_S_2->fill(sphere,cfs.particles().size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale( _h_lam_p    , crossSection()/nanobarn/sumOfWeights());
      scale( _h_lam_pL   , crossSection()/nanobarn/sumOfWeights());
      scale( _h_lam_pTIn , crossSection()/nanobarn/sumOfWeights());
      scale( _h_lam_pTOut, crossSection()/nanobarn/sumOfWeights());
      scale( _h_lam_rap  , crossSection()/nanobarn/sumOfWeights());
      scale( _h_lam_x    , sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      Scatter2DPtr temp;
      book(temp,15+_ih,1,3);
      divide(_p_lam_S_1,_p_lam_S_2,temp);
      if(_ih==0) {
      	scale( _h_xi_p    , crossSection()/nanobarn/sumOfWeights());
      	scale( _h_xi_pL   , crossSection()/nanobarn/sumOfWeights());
      	scale( _h_xi_pTIn , crossSection()/nanobarn/sumOfWeights());
      	scale( _h_xi_pTOut, crossSection()/nanobarn/sumOfWeights());
      	scale( _h_xi_rap  , crossSection()/nanobarn/sumOfWeights());
      	scale( _h_xi_x    , sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_lam_p, _h_lam_pL, _h_lam_pTIn, _h_lam_pTOut, _h_lam_rap, _h_lam_x;
    Profile1DPtr _p_lam_S_1, _p_lam_S_2;
    Histo1DPtr _h_xi_p, _h_xi_pL, _h_xi_pTIn, _h_xi_pTOut, _h_xi_rap, _h_xi_x;
    int _ih;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1989_I266893);


}

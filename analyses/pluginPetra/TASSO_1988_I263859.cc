// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief event shapes at 35 GeV
  class TASSO_1988_I263859 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1988_I263859);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      // Thrust and sphericity
      declare(Beam(), "Beams");
      declare(Thrust(cfs), "Thrust");
      declare(Sphericity(cfs), "Sphericity");

      // Book histograms
      book(_h_sphericity,  1, 1, 1);
      book(_h_aplanarity,  2, 1, 1);
      book(_h_thrust    ,  3, 1, 1);
      book(_h_pTin2     ,  4, 1, 1);
      book(_h_pTout2    ,  5, 1, 1);
      book(_h_ncharged  ,  6, 1, 1);
      book(_h_pTin      ,  7, 1, 1);
      book(_h_pTout     ,  8, 1, 1);
      book(_h_pT        ,  9, 1, 1);
      book(_h_x         , 10, 1, 1);
      book(_h_rap       , 11, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
				   beams.second.p3().mod() ) / 2.0;
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      _h_sphericity->fill(sphericity.sphericity());
      _h_aplanarity->fill(sphericity.aplanarity());
      _h_thrust    ->fill(thrust.thrust());

      double pTin2sum(0.), pTout2sum(0.);
      for (const Particle& p : cfs.particles()) {
        // Get momentum and energy of each particle.
        const Vector3 mom3 = p.p3();
        const double energy = p.E();
        // Scaled momenta.
        const double mom = mom3.mod();
        const double scaledMom = mom/meanBeamMom;
        _h_x->fill(scaledMom);

        const double momS = dot(sphericity.sphericityAxis(), mom3);
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
        const double pT = sqrt(pow(pTinS, 2) + pow(pToutS, 2));
        const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));

	pTin2sum  += sqr(pTinS);
	pTout2sum += sqr(pToutS);

	_h_pTin ->fill(abs(pTinS)/GeV );
	_h_pTout->fill(abs(pToutS)/GeV);
	_h_pT   ->fill(pT/GeV         );
	_h_rap  ->fill(abs(rapidityS) );

	
      }
      unsigned int nCharged = cfs.particles().size();
      _h_pTin2 ->fill(pTin2sum /nCharged);
      _h_pTout2->fill(pTout2sum/nCharged);
      _h_ncharged->fill(nCharged);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sphericity , 1./sumOfWeights());
      scale(_h_aplanarity , 1./sumOfWeights());
      scale(_h_thrust     , 1./sumOfWeights());
      scale(_h_pTin2      , 1./sumOfWeights());
      scale(_h_pTout2     , 1./sumOfWeights());
      scale(_h_ncharged   , 2000./sumOfWeights());
      scale(_h_pTin       , 1./sumOfWeights());
      scale(_h_pTout      , 1./sumOfWeights());
      scale(_h_pT         , 1./sumOfWeights());
      scale(_h_x          , 1./sumOfWeights());
      scale(_h_rap        , 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_sphericity, _h_aplanarity, _h_thrust, _h_pTin2, _h_pTout2, _h_ncharged,
      _h_pTin, _h_pTout, _h_pT, _h_x, _h_rap;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1988_I263859);


}

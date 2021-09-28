// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda and Lambda bar dists
  class DELPHI_1993_I360638 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1993_I360638);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      declare(UnstableParticles(), "UFS");
      declare(Sphericity(cfs), "Sphericity");

      // Book histograms
      book(_h_x       , 1, 1, 1);
      book(_h_rap     , 3, 1, 1);
      book(_h_cos     , 4, 1, 1);
      book(_m_single  , 2, 1, 1);
      book(_m_like    , 5, 1, 1);
      book(_m_opposite, 6, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) vetoEvent;
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      // lambda
      Particles lambda    = ufs.particles(Cuts::pid== PID::LAMBDA);
      Particles lambdabar = ufs.particles(Cuts::pid==-PID::LAMBDA);
      // multiplicities
      _m_single->fill(91.2,(lambda.size()+lambdabar.size()));
      if(lambda.empty()&&lambdabar.empty()) vetoEvent;
      for(const Particle& p : lambda) {
	double xP = 2.*p.p3().mod()/sqrtS();
	_h_x->fill(xP);
      }
      for(const Particle& p : lambdabar) {
	double xP = 2.*p.p3().mod()/sqrtS();
	_h_x->fill(xP);
      }
      if(lambda.size()>=2) {
	unsigned int npair=lambda.size()/2;
	_m_like->fill(91.2,double(npair));
      }
      if(lambdabar.size()>=2) {
	unsigned int npair=lambdabar.size()/2;
	_m_like->fill(91.2,double(npair));
      }
      if(lambda.size()==0 || lambdabar.size()==0)
	return;
      _m_opposite->fill(91.2,double(max(lambda.size(),lambdabar.size())));
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      for(const Particle & p : lambda) {
        const Vector3 momP = p.p3();
        const double  enP  = p.E();
        const double  modP = dot(sphericity.sphericityAxis(), momP);
        const double rapP = 0.5 * std::log((enP + modP) / (enP - modP));
	for(const Particle & pb : lambdabar) {
	  const Vector3 momB = pb.p3();
	  const double  enB  = pb.E();
	  const double  modB = dot(sphericity.sphericityAxis(), momB);
	  const double rapB = 0.5 * std::log((enB + modB) / (enB - modB));
	  _h_rap->fill(abs(rapP-rapB));
	  _h_cos->fill(momP.unit().dot(momB.unit()));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale( _h_x       , 1./sumOfWeights());
      scale( _h_rap     , 1./sumOfWeights());
      scale( _h_cos     , 1./sumOfWeights());
      scale( _m_single  , 1./sumOfWeights());
      scale( _m_like    , 1./sumOfWeights());
      scale( _m_opposite, 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x, _h_rap ,_h_cos;
    Histo1DPtr _m_single, _m_like, _m_opposite;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1993_I360638);


}

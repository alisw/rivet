// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CLEO_2000_I526554 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEO_2000_I526554);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      
      // Book histograms
      book(_h_Ds_star1  , 1, 1, 1);
      book(_h_Ds        , 2, 1, 1);
      book(_h_Ds_star2  , 3, 1, 1);
      book(_h_Ds_primary, 4, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop through unstable FS particles and look for charmed mesons
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");

      const Beam beamproj = apply<Beam>(event, "Beams");
      const ParticlePair& beams = beamproj.beams();
      const FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      LorentzTransform cms_boost;
      if (mom_tot.p3().mod() > 1*MeV)
        cms_boost = LorentzTransform::mkFrameTransformFromBeta(mom_tot.betaVec());
      const double s = sqr(beamproj.sqrtS());
      for (const Particle& p : ufs.particles(Cuts::abspid==431 or Cuts::abspid==433)) {
      	// 3-momentum in CMS frame
      	const double mom = cms_boost.transform(p.momentum()).vector3().mod();
        const int pdgid = p.abspid();
      	double mH2(0.),xp(0.);
      	bool primary = true;
        switch (pdgid) {
      	case 431:
      	  //   MSG_DEBUG("D_s found");
      	  mH2 = sqr(1.96834);
      	  xp = mom/sqrt(s/4.0 - mH2);
      	  _h_Ds->fill(xp);
      	  for(const Particle & mother : p.parents()) {
      	    if(PID::isCharmMeson(mother.pid())) {
      	      primary = false;
      	      break;
      	    }
      	  }
      	  if(primary)
      	    _h_Ds_primary->fill(xp);
      	  break;
        case 433:
          MSG_DEBUG("D_s* found");
          mH2 = sqr(2.1122);
      	  xp = mom/sqrt(s/4.0 - mH2);
      	  _h_Ds_star1->fill(xp);
      	  _h_Ds_star2->fill(xp);
      	  break;
      	default:
      	  break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // BR(D_s->phi pi) x BR(phi->K+K-)
      double br2 = 0.0227;
      // x D_s* -> gamma d_s
      double br1 = br2*.935;
      // cross section factor
      double fact = crossSection()/picobarn/sumOfWeights();
      // normalize the cross sections (bin width)
      scale(_h_Ds_star1  , br1*fact*0.03 );
      scale(_h_Ds        , br2*fact*0.03 );
      scale(_h_Ds_star2  , br2*fact*0.03 );
      scale(_h_Ds_primary, br2*fact*0.03 );

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Ds,_h_Ds_primary,_h_Ds_star1,_h_Ds_star2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CLEO_2000_I526554);


}

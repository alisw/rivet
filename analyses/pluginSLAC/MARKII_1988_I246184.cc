// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  /// @brief Event shapes at 29 GeV
  class MARKII_1988_I246184 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1988_I246184);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      Sphericity sphere(fs);
      declare(sphere, "Sphericity");
      declare(Thrust    (fs), "Thrust"    );
      declare(Hemispheres(sphere), "Hemispheres");
      // histograms
      unsigned int ioff=18;
      for(unsigned int ix=0;ix<3;++ix) {
	book(_histAplanarity [ix]  , 1+ioff*ix, 1, 1);
	book(_histQx         [ix]  , 2+ioff*ix, 1, 1);
	book(_histQ2Q1       [ix]  , 3+ioff*ix, 1, 1);
	book(_histSphericity [ix]  , 4+ioff*ix, 1, 1);
	book(_histThrust     [ix]  , 5+ioff*ix, 1, 1);
	book(_histMinor      [ix]  , 6+ioff*ix, 1, 1);
	book(_histOblateness [ix]  , 7+ioff*ix, 1, 1);
	book(_histMJetBroad  [ix]  , 8+ioff*ix, 1, 1);
	book(_histMJetSlim   [ix]  , 9+ioff*ix, 1, 1);
	book(_histMJetDiff   [ix]  ,10+ioff*ix, 1, 1);
	book(_histScaledMom  [ix]  ,15+ioff*ix, 1, 1);
	book(_histPt2S       [ix]  ,11+ioff*ix, 1, 1);
	book(_histPtS        [ix]  ,12+ioff*ix, 1, 1);
	book(_histPtSIn      [ix]  ,14+ioff*ix, 1, 1);
	book(_histPtSOut     [ix]  ,13+ioff*ix, 1, 1);
	book(_histRapidityS  [ix]  ,16+ioff*ix, 1, 1);
	book(_histTheta      [ix]  ,17+ioff*ix, 1, 1);
	book(_histETheta     [ix]  ,18+ioff*ix, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      for(unsigned int ix=0;ix<3;++ix) {
	_histSphericity[ix]->fill(sphericity.sphericity());
	_histAplanarity[ix]->fill(sphericity.aplanarity());
	_histQx        [ix]->fill((sphericity.lambda1()-sphericity.lambda2())/sqrt(3.));
	_histQ2Q1      [ix]->fill(sphericity.lambda2()-sphericity.lambda3());
      }
      // thrust related
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      for(unsigned int ix=0;ix<3;++ix) {
	_histThrust    [ix]->fill(thrust.thrust());
	_histMinor     [ix]->fill(thrust.thrustMinor());
	_histOblateness[ix]->fill(thrust.oblateness());
      }
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      double mWide = hemi.scaledM2high(), mNarrow = hemi.scaledM2low();
      if(!hemi.massMatchesBroadening()) swap(mWide,mNarrow);
      for(unsigned int ix=0;ix<3;++ix) {
	_histMJetBroad[ix]->fill(mWide);
	_histMJetSlim [ix]->fill(mNarrow);
	_histMJetDiff [ix]->fill(hemi.scaledM2diff());
      }
      // dists w.r.t sphericity axis
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      for (const Particle& p : cfs.particles()) {
        // Get momentum and energy of each particle.
        const Vector3 mom3 = p.p3();
        const double energy = p.E();
        // Scaled momenta.
        const double mom = mom3.mod();
        const double scaledMom = 2.*mom/sqrtS();
	// Get momenta components w.r.t. thrust and sphericity.
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
	double pT2 = sqr(pTinS)+sqr(pToutS);
	double pT  = sqrt(pT2);
        const double momS = dot(sphericity.sphericityAxis(), mom3);
        const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));
	// angle
	double theta = sphericity.sphericityAxis().angle(mom3)/M_PI*180.;
	if(theta>90.) theta=180.-theta;
	// fill histos
	for(unsigned int ix=0;ix<3;++ix) {
	  _histScaledMom[ix]->fill(scaledMom);
	  _histPt2S     [ix]->fill(fabs(pT2/GeV));
	  _histPtS      [ix]->fill(fabs(pT/GeV));
	  _histPtSIn    [ix]->fill(fabs(pTinS/GeV));
	  _histPtSOut   [ix]->fill(fabs(pToutS/GeV));
	  _histRapidityS[ix]->fill(fabs(rapidityS));
	  _histTheta    [ix]->fill(theta);
	}
      }
      // energy flow includes neutral w.r.t sphericity axis
      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles()) {
        // Get momentum and energy of each particle.
        const Vector3 mom3 = p.p3();
        const double energy = p.E();
	// angle
	double theta = sphericity.sphericityAxis().angle(mom3)/M_PI*180.;
	if(theta>90.) theta=180.-theta;
	// fill histos
	for(unsigned int ix=0;ix<3;++ix) {
	  _histETheta   [ix]->fill(theta,energy);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_histAplanarity [ix]  ,1./sumOfWeights());
	scale(_histQx         [ix]  ,1./sumOfWeights());
	scale(_histQ2Q1       [ix]  ,1./sumOfWeights());
	scale(_histSphericity [ix]  ,1./sumOfWeights());
	scale(_histThrust     [ix]  ,1./sumOfWeights());
	scale(_histMinor      [ix]  ,1./sumOfWeights());
	scale(_histOblateness [ix]  ,1./sumOfWeights());
	scale(_histMJetBroad  [ix]  ,1./sumOfWeights());
	scale(_histMJetSlim   [ix]  ,1./sumOfWeights());
	scale(_histMJetDiff   [ix]  ,1./sumOfWeights());
	scale(_histScaledMom  [ix]  ,1./sumOfWeights());
	scale(_histPt2S       [ix]  ,1./sumOfWeights());
	scale(_histPtS        [ix]  ,1./sumOfWeights());
	scale(_histPtSIn      [ix]  ,1./sumOfWeights());
	scale(_histPtSOut     [ix]  ,1./sumOfWeights());
	scale(_histRapidityS  [ix]  ,1./sumOfWeights());
	scale(_histTheta      [ix]  ,1./sumOfWeights());
	scale(_histETheta     [ix]  ,1./sumOfWeights());
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _histAplanarity[3],_histQx[3],_histQ2Q1[3],_histSphericity[3];
    Histo1DPtr _histThrust[3],_histMinor[3],_histOblateness[3];
    Histo1DPtr _histMJetBroad[3],_histMJetSlim[3],_histMJetDiff[3];
    Histo1DPtr _histScaledMom[3],_histPt2S[3],_histPtS[3],_histPtSIn[3],_histPtSOut[3],_histRapidityS[3];
    Histo1DPtr _histTheta[3],_histETheta[3];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1988_I246184);

}

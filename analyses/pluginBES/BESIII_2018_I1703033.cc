// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief omega -> pi+pi-pi0
  class BESIII_2018_I1703033 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1703033);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 223);
      declare(ufs, "UFS");
      DecayedParticles OMEGA(ufs);
      OMEGA.addStable(PID::PI0);
      declare(OMEGA,"OMEGA");
      // histograms
      book(_h_z  ,1,1,1);
      book(_h_phi,1,1,2);
      book(_h_cos,1,1,3);
      book(_h_ppi,1,1,4);
      book(_dalitz,"dalitz",50,-1,1,50,-1.1,0.9);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 211,1}, {-211,1}, { 111,1}};
      static const std::complex<double> ii(0.,1.);
      DecayedParticles OMEGA = apply<DecayedParticles>(event, "OMEGA");
      // loop over particles
      for(unsigned int ix=0;ix<OMEGA.decaying().size();++ix) {
	// pi+ pi- pi0
	if (!OMEGA.modeMatches(ix,3,mode)) continue;
	// apply omega mass cut
	if(abs(OMEGA.decaying()[ix].mass()-0.78265)>0.04) continue;
	// decay products
	const Particle & pip = OMEGA.decayProducts()[ix].at( 211)[0];
	const Particle & pim = OMEGA.decayProducts()[ix].at(-211)[0];
	const Particle & pi0 = OMEGA.decayProducts()[ix].at( 111)[0];
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(OMEGA.decaying()[ix].momentum().betaVec());
	FourMomentum pp = boost.transform(pip.momentum());
	FourMomentum pm = boost.transform(pim.momentum());
	FourMomentum p0 = boost.transform(pi0.momentum());
	// variables from eqn 1 in paper (use version from Eqn 2.1 of 1010.3946)
	double Qc = OMEGA.decaying()[ix].mass()-pi0.mass()-pim.mass()-pip.mass();
	double x = sqrt(3.)*(pm.E()-pp.E())/Qc;
	double y = 3./Qc*(p0.E()-pi0.mass())-1.;
	_dalitz->fill(x,y);
	// plot arg and phase of variable
	std::complex<double> zz = x+ii*y;
	_h_z->fill(norm(zz));
	double phi = arg(zz);
	if(phi<0.) phi+=2.*M_PI;
	_h_phi->fill(phi);
	double ctheta = -pp.p3().unit().dot(p0.p3().unit());
	_h_cos->fill(ctheta);
	_h_ppi->fill(p0.p3().mod());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_z  );
      normalize(_h_phi);
      normalize(_h_cos);
      normalize(_h_ppi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_z,_h_phi,_h_cos,_h_ppi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1703033);

}

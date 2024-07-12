// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi and psi(2S) in Upsilon(1S) decays
  class BELLE_2016_I1454405 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2016_I1454405);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==553), "UFS");
      //histos
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1,1+ix,1+iy);
      book(_nUps,"TMP/nUps");
    }

   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
	if(id == 443 || id==100443) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      for(const Particle & ups : apply<UnstableParticles>(event, "UFS").particles()) {
	_nUps->fill();
	// boost to rest frame
	LorentzTransform boost;
	if (ups.p3().mod() > 1*MeV)
	  boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	// find  psi/psi(2S)
	Particles unstable;
	findDecayProducts(ups,unstable);
	// loop over psi/psi(2S)
	for(const Particle& psi : unstable) {
	  double Pmax = 0.5*(sqr(ups.mass())-sqr(psi.mass()))/ups.mass();
	  const FourMomentum p2 = boost.transform(psi.momentum());
	  const double xP = p2.p3().mod()/Pmax;
	  if(psi.pid()==443) {
	    _h[0][0]->fill(xP);
	    _h[1][0]->fill(xP);
	  }
	  else {
	    _h[0][1]->fill(xP);
	    _h[1][1]->fill(xP);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  scale(_h[ix][iy],1e4/ *_nUps);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    CounterPtr _nUps;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2016_I1454405);

}

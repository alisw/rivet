// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c1 in Upsilon(1,2S) decays
  class BELLE_2017_I1499706 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2017_I1499706);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==553 || Cuts::pid==100553), "UFS");
      //histos
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1,1+ix,1+iy);
	book(_c[ix],"TMP/nUps"+toString(ix+1));
      }
    }

   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
	if(id == 20443) {
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
	unsigned int iloc = ups.pid()/100000;
	_c[iloc]->fill();
	// boost to rest frame
	LorentzTransform boost;
	if (ups.p3().mod() > 1*MeV)
	  boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	// find ch_c1
	Particles unstable;
	findDecayProducts(ups,unstable);
	// loop over chi_c1
	for(const Particle& chi : unstable) {
	  double Pmax = 0.5*(sqr(ups.mass())-sqr(chi.mass()))/ups.mass();
	  const FourMomentum p2 = boost.transform(chi.momentum());
	  const double xP = p2.p3().mod()/Pmax;
	  _h[0][iloc]->fill(xP);
	  _h[1][iloc]->fill(xP);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  scale(_h[ix][iy],1e4/ *_c[iy]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2017_I1499706);

}

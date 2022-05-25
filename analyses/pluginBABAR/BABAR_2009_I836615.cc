// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D*+/- in Upsilon(1S) decays
  class BABAR_2009_I836615 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I836615);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // book the histogram
      book(_hist,1,1,1);
    }

   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
	if(id == 413) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // maximum momentum of D* calculated using PDG 2008 masses (as in paper)
      static const double Pmax = 4.28172;
      // Loop over the upsilons
      // Find the Upsilons among the unstables
      for(const Particle & ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==553)) {
	// boost to rest frame
	LorentzTransform boost;
	if (ups.p3().mod() > 1*MeV)
	  boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	// find D*
	Particles unstable;
	findDecayProducts(ups,unstable);
	// loop over D*
	for(const Particle& p : unstable) {
	  const FourMomentum p2 = boost.transform(p.momentum());
	  const double xP = p2.p3().mod()/Pmax;
	  _hist->fill(xP);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_hist);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _hist;
    //@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I836615);

}

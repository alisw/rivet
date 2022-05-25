// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief D_1 D*2 decay angle
  class CLEO_1990_I281039 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1990_I281039);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_D1, 1, 1, 1);
      book(_h_D2, 1, 1, 2);

    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles & dstar, Particles & pi,unsigned int & ncount) {
      for(const Particle & p: mother.children()) {
	if(p.abspid()==413)
	  dstar.push_back(p);
	else if(p.abspid()==211)
	  pi.push_back(p);
	ncount +=1;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==425||
									       Cuts::abspid==10423)) {
	// decay products
	Particles dstar,pi;
	unsigned int ncount=0;
	findDecayProducts(p,dstar,pi,ncount);
	if(ncount!=2 || dstar.size()!=1 || pi.size()!=1) continue;
	if(dstar[0].pid()/p.pid()<0) continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 axis = boost.transform(dstar[0].momentum()).p3().unit();
	double cosL  = axis.dot(p.momentum().p3().unit());
	// decay angles
	if(p.abspid()==425)
	  _h_D2->fill(cosL);
	else
	  _h_D1->fill(cosL);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_D1);
      normalize(_h_D2);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_D1,_h_D2;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CLEO_1990_I281039);


}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_1 and D_2 spectra and decay distributions
  class ARGUS_1989_I280943 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1989_I280943);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_D1_rate    ,1,1,1);
      book(_h_D2_rate    ,1,1,2);
      book(_h_D1_x       ,4,1,1);
      book(_h_D2_x       ,4,1,2);
      book(_h_D1_alpha   ,3,1,1);
      book(_h_D2_alpha   ,3,1,2);
    }
    
    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles & dstar, Particles & d0, Particles & pi,unsigned int & ncount) {
      for(const Particle & p: mother.children()) {
	if(p.abspid()==413)
	  dstar.push_back(p);
	else if(p.abspid()==421)
	  d0.push_back(p);
	else if(p.abspid()==211)
	  pi.push_back(p);
	ncount +=1;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==425 || Cuts::abspid==10423)) {
	const double xp = 2.*p.p3().mod()/sqrtS();
	// spectra
	if(p.abspid()==425)
	  _h_D2_x->fill(xp);
	else
	  _h_D1_x->fill(xp);
	// decay products
	// first od D_1,D_2
	Particles dstar,d0,pi;
	unsigned int ncount=0;
	findDecayProducts(p,dstar,d0, pi,ncount);
	if(ncount!=2 || dstar.size()!=1 || pi.size()!=1 || d0.size()!=0 ) continue;
	if(dstar[0].pid()/p.pid()<0) continue;
	if(p.abspid()==425)
	  _h_D2_rate->fill(10.);
	else
	  _h_D1_rate->fill(10.);
	Particle p2 = dstar[0];
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p2.momentum().betaVec());
	Vector3 d1 = boost.transform(pi[0].momentum()).p3().unit();
	// then of D*
	ncount=0;
	dstar.clear();
	d0.clear();
	pi.clear();
	findDecayProducts(p2,dstar,d0, pi,ncount);
	if(ncount!=2 || dstar.size()!=0 || pi.size()!=1 || d0.size()!=1 ) continue;
	if(pi[0].pid()/p2.pid()<0) continue;
	Vector3 d2 = boost.transform(pi[0].momentum()).p3().unit();
	double cosAlpha  = abs(d1.dot(d2));
	// decay angles
	if(p.abspid()==425)
	  _h_D2_alpha->fill(cosAlpha);
	else
	  _h_D1_alpha->fill(cosAlpha);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_D1_x);
      normalize(_h_D2_x);
      normalize(_h_D1_alpha);
      normalize(_h_D2_alpha);
      scale(_h_D1_rate,crossSection()/picobarn/sumOfWeights());
      scale(_h_D2_rate,crossSection()/picobarn/sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_D1_rate, _h_D2_rate;
    Histo1DPtr _h_D1_x, _h_D2_x;
    Histo1DPtr _h_D1_alpha, _h_D2_alpha;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ARGUS_1989_I280943);


}

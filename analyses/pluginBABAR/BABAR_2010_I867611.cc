// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  D_1 and D_2 decay distributions
  class BABAR_2010_I867611 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2010_I867611);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_D1_ctheta, 1,1,1);
      book(_h_D2_ctheta, 1,1,2);
      book(_h_D02S_ctheta, 1,1,3);
      book(_h_DStar02S_ctheta, 1,1,4);
      book(_h_D3_ctheta, 1,1,5);
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
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==425 || Cuts::abspid==10423 ||
									       Cuts::abspid==100421 || Cuts::abspid==100423 ||
									       Cuts::abspid==427)) {
	// decay products
	Particles dstar,d0,pi;
	unsigned int ncount=0;
	findDecayProducts(p,dstar,d0, pi,ncount);
	if(ncount!=2 || dstar.size()!=1 || pi.size()!=1 || d0.size()!=0 ) continue;
	if(dstar[0].pid()/p.pid()<0) continue;
	Particle p2 = dstar[0];
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p2.momentum().betaVec());
	Vector3 d1 = boost.transform(pi[0].momentum()).p3().unit();
	ncount=0;
	dstar.clear();
	d0.clear();
	pi.clear();
	findDecayProducts(p2,dstar,d0, pi,ncount);
	if(ncount!=2 || dstar.size()!=0 || pi.size()!=1 || d0.size()!=1 ) continue;
	if(pi[0].pid()/p2.pid()<0) continue;
	Vector3 d2 = boost.transform(pi[0].momentum()).p3().unit();
	double cTheta  = d1.dot(d2);
	// decay angles
	if(p.abspid()==425)
	  _h_D2_ctheta->fill(cTheta);
	else if(p.abspid()==10423) 
	  _h_D1_ctheta->fill(cTheta);
	else if(p.abspid()==100421) 
	  _h_D02S_ctheta->fill(cTheta);
	else if(p.abspid()==100423) 
	  _h_DStar02S_ctheta->fill(cTheta);
	else if(p.abspid()==427) 
	  _h_D3_ctheta->fill(cTheta);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_D1_ctheta);
      normalize(_h_D2_ctheta);
      normalize(_h_D02S_ctheta);
      normalize(_h_DStar02S_ctheta);
      normalize(_h_D3_ctheta);
    }

    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_D1_ctheta,_h_D2_ctheta,_h_D02S_ctheta;
    Histo1DPtr _h_DStar02S_ctheta,_h_D3_ctheta;
    //@}

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2010_I867611);


}

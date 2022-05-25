// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D**_s decays
  class BABAR_2009_I827985 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I827985);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_DStar_ctheta, 1,1,1);
      book(_h_D3_ctheta[0], 1,1,2);
      book(_h_D3_ctheta[1], 1,1,3);
      book(_h_D3_ctheta[2], 1,1,4);
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles & dstar, Particles & d0, Particles & K0, Particles & pi, unsigned int & ncount) {
      for(const Particle & p: mother.children()) {
	if(p.abspid()==413)
	  dstar.push_back(p);
	else if(p.abspid()==421)
	  d0.push_back(p);
	else if(p.abspid()==130 || p.abspid()==130 || p.abspid()==311)
	  K0.push_back(p);
	else if(p.abspid()==211)
	  pi.push_back(p);
	ncount +=1;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==100433 || Cuts::abspid==437 ||
									       Cuts::abspid==30433)) {
	// decay products
	Particles dstar,d0,K0,pi;
	unsigned int ncount=0;
	findDecayProducts(p, dstar, d0, K0, pi, ncount);
	if(ncount!=2 || dstar.size()!=1 || K0.size()!=1 ) continue;
	if(dstar[0].pid()/p.pid()<0) continue;
	Particle p2 = dstar[0];
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p2.momentum().betaVec());
	Vector3 d1 = boost.transform(K0[0].momentum()).p3().unit();
	ncount=0;
	dstar.clear();
	d0.clear();
	pi.clear();
	findDecayProducts(p2, dstar, d0, K0, pi, ncount);
	if(ncount!=2 || pi.size()!=1 || d0.size()!=1 ) continue;
	if(pi[0].pid()/p2.pid()<0) continue;
	Vector3 d2 = boost.transform(pi[0].momentum()).p3().unit();
	double cTheta  = d1.dot(d2);
	// decay angles
	if(p.abspid()==100433)
	  _h_DStar_ctheta->fill(cTheta);
	else if(p.abspid()==30433) {
	  _h_D3_ctheta[0]->fill(cTheta);
	  _h_D3_ctheta[2]->fill(cTheta);
	}
	else if(p.abspid()==437) {
	  _h_D3_ctheta[1]->fill(cTheta);
	  _h_D3_ctheta[2]->fill(cTheta);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_DStar_ctheta);
      normalize(_h_D3_ctheta[0]);
      normalize(_h_D3_ctheta[1]);
      normalize(_h_D3_ctheta[2]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_DStar_ctheta,_h_D3_ctheta[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I827985);

}

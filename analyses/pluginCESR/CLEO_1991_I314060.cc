// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief D*+ polarization
  class CLEO_1991_I314060 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1991_I314060);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      
      {Histo1DPtr temp; _h_ctheta.add(0.25,0.45,book(temp,2,1,1));}
      {Histo1DPtr temp; _h_ctheta.add(0.45,0.55,book(temp,2,1,2));}
      {Histo1DPtr temp; _h_ctheta.add(0.55,0.65,book(temp,2,1,3));}
      {Histo1DPtr temp; _h_ctheta.add(0.65,0.75,book(temp,2,1,4));}
      {Histo1DPtr temp; _h_ctheta.add(0.75,0.85,book(temp,2,1,5));}
      {Histo1DPtr temp; _h_ctheta.add(0.85,1.00,book(temp,2,1,6));}

    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles & d0, Particles & pi,unsigned int & ncount) {
      for(const Particle & p: mother.children()) {
	if(p.abspid()==421)
	  d0.push_back(p);
	else if(p.abspid()==211)
	  pi.push_back(p);
	ncount +=1;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // loop over D*+ mesons
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p: ufs.particles(Cuts::abspid==413)) {
	// calc x+
	double x = (p.momentum().E()+p.momentum().z())/(meanBeamMom + sqrt(sqr(meanBeamMom)+p.mass2()));
	// checck decay products
	Particles d0,pi;
	unsigned int ncount=0;
	findDecayProducts(p,d0, pi,ncount);
	if(ncount!=2 || pi.size()!=1 || d0.size()!=1 ) continue;
	if(pi[0].pid()/p.pid()<0) continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 d1 = boost.transform(pi[0].momentum()).p3().unit();
	double ctheta  = d1.dot(p.momentum().p3().unit());
	_h_ctheta.fill(x,ctheta);
      }
    }

    pair<double,pair<double,double> > calcAlpha(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,make_pair(0.,0.));
      double sum1(0.),sum2(0.),sum3(0.),sum4(0.),sum5(0.);
      for (auto bin : hist->bins() ) {
       	double Oi = bin.area();
	if(Oi==0.) continue;
	double a =  1.5*(bin.xMax() - bin.xMin());
	double b = 0.5*(pow(bin.xMax(),3) - pow(bin.xMin(),3));
       	double Ei = bin.areaErr();
	sum1 +=   a*Oi/sqr(Ei);
	sum2 +=   b*Oi/sqr(Ei);
	sum3 += sqr(a)/sqr(Ei);
	sum4 += sqr(b)/sqr(Ei);
	sum5 +=    a*b/sqr(Ei);
      }
      // calculate alpha
      double alpha = (-3*sum1 + 9*sum2 + sum3 - 3*sum5)/(sum1 - 3*sum2 + 3*sum4 - sum5);
      // and error
      double cc = -pow((sum3 + 9*sum4 - 6*sum5),3);
      double bb = -2*sqr(sum3 + 9*sum4 - 6*sum5)*(sum1 - 3*sum2 + 3*sum4 - sum5);
      double aa =  sqr(sum1 - 3*sum2 + 3*sum4 - sum5)*(-sum3 - 9*sum4 + sqr(sum1 - 3*sum2 + 3*sum4 - sum5) + 6*sum5);      
      double dis = sqr(bb)-4.*aa*cc;
      if(dis>0.) {
	dis = sqrt(dis);
	return make_pair(alpha,make_pair(0.5*(-bb+dis)/aa,-0.5*(-bb-dis)/aa));
      }
      else {
	return make_pair(alpha,make_pair(0.,0.));
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      vector<double> x   = {0.35,0.5 ,0.6 ,0.7 ,0.8 ,0.925};
      vector<double> wid = {0.10,0.05,0.05,0.05,0.05,0.075};
      Scatter2DPtr h_alpha;
      book(h_alpha,1,1,1);
      Scatter2DPtr h_rho  ;
      book(h_rho  ,1,1,2);
      Scatter2DPtr h_eta  ;
      book(h_eta  ,1,1,3);
      for(unsigned int ix=0;ix<_h_ctheta.histos().size();++ix) {
	// normalize
     	normalize(_h_ctheta.histos()[ix]);
	// alpha
    	pair<double,pair<double,double> > alpha = calcAlpha(_h_ctheta.histos()[ix]);
	h_alpha->addPoint(x[ix], alpha.first, make_pair(wid[ix],wid[ix]),alpha.second);
	// rho
	double rho = (1.+alpha.first)/(3.+alpha.first);
	pair<double,double> rho_error;
	rho_error.first  = 2.*alpha.second.first /sqr(3.+alpha.first);
	rho_error.second = 2.*alpha.second.second/sqr(3.+alpha.first);
	h_rho->addPoint(x[ix], rho, make_pair(wid[ix],wid[ix]),rho_error);
	// eta
	double eta = alpha.first/(3.+alpha.first);
	pair<double,double> eta_error;
	eta_error.first  = 3.*alpha.second.first /sqr(3.+alpha.first);
	eta_error.second = 3.*alpha.second.second/sqr(3.+alpha.first);
	h_eta->addPoint(x[ix], eta, make_pair(wid[ix],wid[ix]),eta_error);
      }
    }

    //@}


    /// @name Histograms
    //@{
    BinnedHistogram _h_ctheta;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CLEO_1991_I314060);


}

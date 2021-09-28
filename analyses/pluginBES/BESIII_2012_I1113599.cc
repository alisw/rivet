// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi to p pbar n nbar
  class BESIII_2012_I1113599 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2012_I1113599);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      declare(FinalState(), "FS");

      // Book histograms
      book(_h_proton , "ctheta_p",20,-1.,1.);
      book(_h_neutron, "ctheta_n",20,-1.,1.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // types of final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      Particle outgoing;
      for (const Particle& p :  fs.particles()) {
	nCount[p.pid()] += 1;
	if(p.pid()==2212 || p.pid()==2112)
	  outgoing = p;
	++ntotal;
      }
      if(ntotal==2) {
	if(nCount[2212]==1 && nCount[-2212]==1) {
	  _h_proton->fill(outgoing.momentum().p3().unit().dot(axis));
	}
	else if(nCount[2112]==1 && nCount[-2112]==1) {
	  _h_neutron->fill(outgoing.momentum().p3().unit().dot(axis));
	}
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
      // proton
      normalize(_h_proton );
      pair<double,pair<double,double> > alpha = calcAlpha(_h_proton);
      Scatter2DPtr _h_alpha_proton;
      book(_h_alpha_proton,1,1,1);
      _h_alpha_proton->addPoint(0.5, alpha.first, make_pair(0.5,0.5),
				make_pair(alpha.second.first,alpha.second.second) );
      // neutron
      normalize(_h_neutron);
      alpha = calcAlpha(_h_neutron);
      Scatter2DPtr _h_alpha_neutron;
      book(_h_alpha_neutron, 1,1,2);
      _h_alpha_neutron->addPoint(0.5, alpha.first, make_pair(0.5,0.5),
				 make_pair(alpha.second.first,alpha.second.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_proton,_h_neutron;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESIII_2012_I1113599);


}

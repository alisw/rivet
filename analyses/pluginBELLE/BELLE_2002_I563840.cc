// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"


namespace Rivet {


  /// @brief charmonium production at 10.6 GeV
  class BELLE_2002_I563840 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2002_I563840);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // book the histograms
      // spectra
      book(_h_Jpsi,3,1,1);
      book(_h_feed,3,1,2);
      book(_h_Psi2,3,2,1);
      // cross sections
      book(_h_sig_JPsi_all ,1,1,1);
      book(_h_sig_Jpsi_high,1,2,1);
      book(_h_sig_Jpsi_feed,1,2,2);
      book(_h_sig_Psi2_high,1,2,3);
      // angular distributions
      vector<double> bins = {2.,2.6,3.4,4.9};
      for(unsigned int ix=0;ix<3;++ix) {
      	Histo1DPtr temp;
      	// book cos theta*
      	if(ix<=1) {
      	  std::ostringstream title;
      	  title << "/TMP/cThetaStar_" << ix+1;
      	  book(temp,title.str(),5,-1.,1.);
      	}
      	else {
      	  book(temp,4,1,2);
      	}
      	_h_cThetaStar.add(bins[ix],bins[ix]+1,temp);
      	// book cos theta_H
      	if(ix<=1) {
      	  std::ostringstream title;
      	  title << "/TMP/cThetaH_" << ix+1;
      	  book(temp,title.str(),5,-1.,1.);
      	}
      	else {
      	  book(temp,4,2,2);
      	}
      	_h_cThetaH.add(bins[ix],bins[ix]+1,temp);
      }
      book(_h_cS_low,4,1,1);
      book(_h_cS_high,"/TMP/cS_high",5,-1.,1.);
      book(_h_cH_low,4,2,1);
      book(_h_cH_high,"/TMP/cH_high",5,-1.,1.);
    }

    void findLeptons(const Particle & mother,
		     unsigned int & nstable,
		     Particles& lp, Particles& lm) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == 11 || id == 13 ) {
	  lm.push_back(p);
	  ++nstable;
	}
       	else if (id == -11 || id==-13) {
       	  lp.push_back(p);
       	  ++nstable;
       	}
	else if ( !p.children().empty() ) {
	  findLeptons(p,nstable,lp,lm);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & p : apply<UnstableParticles>("UFS",event).particles(Cuts::pid==443 or Cuts::pid==100443 )) {
      	LorentzTransform boost = cmsTransform( beams() );
      	// check if prompt (i.e. not from B decay)
      	if(p.fromBottom()) continue;
      	bool feedDown = false;
	FourMomentum mom = boost.transform(p.momentum());
	double pStar = mom.p3().mod();
      	if(p.pid()==443) {
      	  Particle parent=p;
      	  while(!parent.parents().empty()) {
      	    parent=parent.parents()[0];
      	    if(p.pid()==parent.pid()) continue;
      	    if((parent.abspid()%1000)/10==44) {
      	      feedDown=true;
      	      break;
      	    }
      	  }
	  _h_Jpsi->fill(pStar);
	  _h_sig_JPsi_all->fill(10.6);
	  if(pStar>2.) {
	    _h_sig_Jpsi_high->fill(10.6);
	    if(feedDown) {
	      _h_feed->fill(pStar);
	      _h_sig_Jpsi_feed->fill(10.6);
	    }
	    double cThetaS = cos(mom.p3().polarAngle());
	    _h_cThetaStar.fill(pStar,cThetaS);
	    if(pStar<3.4) _h_cS_low ->fill(cThetaS);
	    else          _h_cS_high->fill(cThetaS);
	    // leptons from J/psi decay
	    unsigned int nstable = 0;
	    Particles lp, lm;
	    findLeptons(p,nstable,lp,lm);
	    if(nstable==2&&lp.size()==1&&lm.size()==1) {
	      FourMomentum pl = boost.transform(lp[0].momentum());
	      LorentzTransform b2 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	      pl = b2.transform(pl);
	      double cThetaH = pl.p3().unit().dot(p.p3().unit());
	      _h_cThetaH.fill(pStar,cThetaH);
	      if(pStar<3.4) _h_cH_low ->fill(cThetaH);
	      else          _h_cH_high->fill(cThetaH);
	    }
	  }
	}
	else {
	  _h_Psi2->fill(pStar);
	  if(pStar>2.) _h_sig_Psi2_high->fill(10.6);
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
      // spectra
      normalize(_h_Jpsi,1.,false);
      normalize(_h_feed,1.,false);
      normalize(_h_Psi2,1.,false);
      // cross sections
      double fact = 1./sumOfWeights()*crossSection()/picobarn;
      scale(_h_sig_JPsi_all ,fact);
      scale(_h_sig_Jpsi_high,fact);
      scale(_h_sig_Jpsi_feed,fact);
      scale(_h_sig_Psi2_high,fact);
      // angular dists and parameters from them
      vector<double> bins = {2.,2.6,3.4,4.9};
      Scatter2DPtr _h_A;
      book(_h_A    ,2,1,1);
      Scatter2DPtr _h_alpha;
      book(_h_alpha,2,1,2);
      for(unsigned int ix=0;ix<3;++ix) {
      	normalize(_h_cThetaStar.histos()[ix]);
      	pair<double,pair<double,double> > alpha = calcAlpha(_h_cThetaStar.histos()[ix]);
      	double centre = 0.5*(bins[ix+1]+bins[ix]);
      	double width  = 0.5*(bins[ix+1]-bins[ix]);
      	_h_A->addPoint(centre, alpha.first, make_pair(width,width),
      		       make_pair(alpha.second.first,alpha.second.second) );
      	normalize(_h_cThetaH.histos()[ix]);
      	alpha = calcAlpha(_h_cThetaH.histos()[ix]);
      	_h_alpha->addPoint(centre, alpha.first, make_pair(width,width),
      			   make_pair(alpha.second.first,alpha.second.second) );
      }
      double centre1 = 0.5*(bins[2]+bins[0]);
      double width1  = 0.5*(bins[2]-bins[0]);
      double centre2 = 0.5*(bins[3]+bins[0]);
      double width2  = 0.5*(bins[3]-bins[0]);
      normalize(_h_cS_low);
      pair<double,pair<double,double> > alpha = calcAlpha(_h_cS_low);
      book(_h_A    ,2,2,1);
      _h_A->addPoint(centre1, alpha.first, make_pair(width1,width1),
      		     make_pair(alpha.second.first,alpha.second.second) );
      normalize(_h_cS_high);
      alpha = calcAlpha(_h_cS_high);
      book(_h_A    ,2,3,1);
      _h_A->addPoint(centre2, alpha.first, make_pair(width2,width2),
      		     make_pair(alpha.second.first,alpha.second.second) );
      normalize(_h_cH_low);
      alpha = calcAlpha(_h_cH_low);
      book(_h_alpha,2,2,2);
      _h_alpha->addPoint(centre1, alpha.first, make_pair(width1,width1),
      		     make_pair(alpha.second.first,alpha.second.second) );
      normalize(_h_cH_high);
      alpha = calcAlpha(_h_cH_high);
      book(_h_alpha,2,3,2);
      _h_alpha->addPoint(centre2, alpha.first, make_pair(width2,width2),
      			 make_pair(alpha.second.first,alpha.second.second) );
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Jpsi,_h_Psi2,_h_feed;
    Histo1DPtr _h_sig_JPsi_all,_h_sig_Jpsi_high,_h_sig_Jpsi_feed,_h_sig_Psi2_high;
    Histo1DPtr _h_cS_low,_h_cS_high,_h_cH_low,_h_cH_high;
    BinnedHistogram _h_cThetaStar,_h_cThetaH;


    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2002_I563840);

}

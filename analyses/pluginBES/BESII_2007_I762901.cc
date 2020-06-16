// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESII_2007_I762901 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESII_2007_I762901);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=1;ix<11;++ix) {
	stringstream ss;
	ss << "TMP/n" << ix;
	book(_nMeson[ix], ss.str());
      }

    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==6) {
	if(nCount[211]==2 && nCount[-211]==2 && nCount[321]==1 && nCount[-321]==1)
	  _nMeson[1]->fill();
	else if(nCount[211]==1 && nCount[-211]==1 && nCount[321]==2 && nCount[-321]==2)
	  _nMeson[2]->fill();
	else if(nCount[211]==2 && nCount[-211]==2 && nCount[2212]==1 && nCount[-2212]==1)
	  _nMeson[3]->fill();
	else if(nCount[211]==3 && nCount[-211]==3)
	  _nMeson[4]->fill();
      }
      else if(ntotal==7) {
	if(nCount[211]==2 && nCount[-211]==2 && nCount[321]==1 && nCount[-321]==1 && nCount[111]==1)
	  _nMeson[5]->fill();
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
       	if(p.children().empty()) continue;
       	if(p.pid()!=113 && abs(p.pid())!=313) continue;
     	map<long,int> nRes = nCount;
     	int ncount = ntotal;
     	findChildren(p,nRes,ncount);
      	if(p.pid()==113) {
       	  if(ncount!=2) continue;
      	  unsigned int nPi(0),nK(0),nProt(0);
       	  bool matched = true;
       	  for(auto const & val : nRes) {
       	    if(abs(val.first)==211 && val.second==1) {
       	      nPi+=1;
       	    }
       	    else if(abs(val.first)==321 && val.second==1) {
       	      nK+=1;
       	    }
       	    else if(abs(val.first)==2212 && val.second==1) {
       	      nProt+=1;
       	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
       	  }
       	  if(matched) {
       	    if(nProt==2 && nPi==0 && nK==0)
       	      _nMeson[9]->fill();
       	    else if(nProt==0 && nPi==2 && nK==0)
       	      _nMeson[7]->fill();
       	    else if(nProt==0 && nPi==0 && nK==2)
       	      _nMeson[8]->fill();
       	  }
       	}
     	else if(abs(p.pid())==313) {
       	  if(ncount!=2) continue;
      	  unsigned int npi(0),nK(0);
      	  bool matched = true;
      	  int ipi = p.pid()==313 ?  211 : -211;
       	  int iK  = p.pid()==313 ? -321 :  321;
      	  for(auto const & val : nRes) {
      	    if(abs(val.first)== ipi && val.second==1) {
      	      npi+=1;
      	    }
      	    else if(abs(val.first)==iK && val.second==1) {
      	      nK+=1;
      	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
       	  }
       	  if(matched) {
       	    if(npi==1&&nK==1)
       	      _nMeson[10]->fill();
       	  }
       	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<11;++ix) {
	if(ix==6) continue;
	double sigma = _nMeson[ix]->val();
	double error = _nMeson[ix]->err();
    	sigma *= crossSection()/ sumOfWeights() /picobarn;
    	error *= crossSection()/ sumOfWeights() /picobarn;
	Scatter2D temphisto(refData(1, 1, ix));
	Scatter2DPtr  mult;
        book(mult, 1, 1, ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, sigma, ex, make_pair(error,error));
	  }
	  else {
	    mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _nMeson[11];
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESII_2007_I762901);


}

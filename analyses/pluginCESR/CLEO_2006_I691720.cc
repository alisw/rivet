// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CLEO_2006_I691720 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2006_I691720);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=1;ix<18;++ix) {
	stringstream ss;
	ss << "TMP/n" << ix;
	book(_nMeson[ix], ss.str());
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      if(ntotal==3) {
	if(nCount[211]==1 && nCount[-211]==1 && nCount[111]==1)
	  _nMeson[1]->fill();
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
       	if(p1.children().empty()) continue;
       	if(p1.pid()!=113 && abs(p1.pid())!=313) continue;
     	map<long,int> nRes = nCount;
     	int ncount = ntotal;
     	findChildren(p1,nRes,ncount);
      	if(p1.pid()==113 || p1.pid()==223 || p1.pid()==10113) {
	  if(ncount!=1) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==111 && val.second!=1) {
	      matched = false;
	      break;
	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
	  }
       	  if(matched) {
	    if(p1.pid()==113) {
	      _nMeson[2]->fill();
	      _nMeson[3]->fill();
	    }
	    else if(p1.pid()==223) {
	      _nMeson[5]->fill();
	    }
	    else if(p1.pid()==10113) {
	      _nMeson[15]->fill();
	    }
	  }
	}
      	else if(p1.pid()==abs(113) || abs(p1.pid())==10113) {
	  if(ncount!=1) continue;
	  bool matched = true;
	  int ipi = p1.pid()>0 ? -211 : 211;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==ipi && val.second!=1) {
	      matched = false;
	      break;
	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
	  }
       	  if(matched) {
	    if(p1.pid()==abs(113)) {
	      _nMeson[2]->fill();
	      _nMeson[4]->fill();
	    }
	    else {
	      _nMeson[15]->fill();
	      _nMeson[17]->fill();
	    }
	  }
	}
      	else if(p1.pid()==abs(323)) {
	  if(ncount!=1) continue;
	  bool matched = true;
	  int iK = p1.pid()==323 ? -321 : 321;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==iK && val.second!=1) {
	      matched = false;
	      break;
	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
	  }
       	  if(matched) {
	    _nMeson[14]->fill();
	  }
	}
	// second unstable particle
	if(abs(p1.pid())!=313 && p1.pid()!=113 && p1.pid()!=223 &&
	   p1.pid()!=333 && p1.pid()!=221 && p1.pid()!=331)
	  continue;
	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
	  const Particle& p2 = ufs.particles()[iy];
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount!=0) continue;
	  bool matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(!matched) continue;
	  if( (p1.pid()==113 && p2.pid()==221) ||
	      (p2.pid()==113 && p1.pid()==221) )
	    _nMeson[7]->fill();
	  else if( (p1.pid()==223 && p2.pid()==221) ||
		   (p2.pid()==223 && p1.pid()==221) )
	    _nMeson[8]->fill();
	  else if( (p1.pid()==333 && p2.pid()==221) ||
		   (p2.pid()==333 && p1.pid()==221) )
	    _nMeson[9]->fill();
	  else  if( (p1.pid()==113 && p2.pid()==331) ||
		    (p2.pid()==113 && p1.pid()==331) )
	    _nMeson[10]->fill();
	  else  if( (p1.pid()==223 && p2.pid()==331) ||
		    (p2.pid()==223 && p1.pid()==331) )
	    _nMeson[11]->fill();
	  else  if( (p1.pid()==333 && p2.pid()==331) ||
		    (p2.pid()==333 && p1.pid()==331) )
	    _nMeson[12]->fill();
	  else  if( (p1.pid()==313 && p2.pid()==-313) ||
		    (p2.pid()==313 && p1.pid()==-313) )
	    _nMeson[13]->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<18;++ix) {
	if(ix==6||ix==16) continue;
	double sigma = _nMeson[ix]->val();
	double error = _nMeson[ix]->err();
    	sigma *= crossSection()/ sumOfWeights() /picobarn;
    	error *= crossSection()/ sumOfWeights() /picobarn;
	Scatter2D temphisto(refData(1, 1, ix));
	Scatter2DPtr mult;
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
    CounterPtr _nMeson[18];
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CLEO_2006_I691720);


}

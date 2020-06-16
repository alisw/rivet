// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CLEOC_2008_I777917 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOC_2008_I777917);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      book(_c_D0D0, "/TMP/sigma_D0D0");
      book(_c_DpDm, "/TMP/sigma_DpDm");
      book(_c_DsDs, "/TMP/sigma_DsDs");
      book(_c_D0D0S, "/TMP/sigma_D0D0S");
      book(_c_DpDmS, "/TMP/sigma_DpDmS");
      book(_c_DsDsS, "/TMP/sigma_DsDsS");
      book(_c_D0SD0S, "/TMP/sigma_D0SD0S");
      book(_c_DpSDmS, "/TMP/sigma_DpSDmS");
      book(_c_DsSDsS, "/TMP/sigma_DsSDsS");
      book(_c_DD, "/TMP/sigma_DD");
      book(_c_DDX, "/TMP/sigma_DDX");
      book(_c_DSDpi, "/TMP/sigma_DSDpi");
      book(_c_DSDSpi, "/TMP/sigma_DSDSpi");
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill();
      // everything else
      else
	_c_hadrons->fill();
      // identified final state with D mesons
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
	int id1 = abs(p1.pid());
	if(id1 != 411 && id1 != 413 && id1 != 421 && id1 != 423 &&
	   id1 != 431 && id1 != 433)
	  continue;
	// check fs
	bool fs = true;
	for (const Particle & child : p1.children()) {
	  if(child.pid()==p1.pid()) {
	    fs = false;
	    break;
	  }
	}
	if(!fs) continue;
	// find the children
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	bool matched=false;
	int sign = p1.pid()/id1;
	// loop over the other fs particles
	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
	  const Particle& p2 = ufs.particles()[iy];
	  fs = true;
	  for (const Particle & child : p2.children()) {
	    if(child.pid()==p2.pid()) {
	      fs = false;
	      break;
	    }
	  }
	  if(!fs) continue;
	  if(p2.pid()/abs(p2.pid())==sign) continue;
	  int id2 = abs(p2.pid());
	  if(id2 != 411 && id2 != 413 && id2 != 421 && id2 != 423 &&
	     id2 != 431 && id2 != 433)
	    continue;
	  if(!p2.parents().empty() && p2.parents()[0].pid()==p1.pid())
	    continue;
	  if((id1==411 || id1==421 || id1==431) && (id2==411 || id2==421 || id2==431 ))
	    _c_DDX->fill();
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2==0) {
	    matched=true;
	    for(auto const & val : nRes2) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      if(id1==411 && id2==411) {
		_c_DpDm->fill();
		_c_DD  ->fill();
	      }
	      else if(id1==421&& id2==421) {
		_c_D0D0->fill();
		_c_DD  ->fill();
	      }
	      else if(id1==431&& id2==431) {
		_c_DsDs->fill();
	      }
	      else if(id1==413 && id2==413) {
		_c_DpSDmS->fill();
	      }
	      else if(id1==423&& id2==423) {
		_c_D0SD0S->fill();
	      }
	      else if(id1==433&& id2==433) {
		_c_DsSDsS->fill();
	      }
	      else if((id1==421 && id2==423) ||
		      (id1==423 && id2==421)) {
		_c_D0D0S->fill();
	      }
	      else if((id1==411 && id2==413) ||
		      (id1==413 && id2==411)) {
		_c_DpDmS->fill();
	      }
	      else if((id1==431 && id2==433) ||
		      (id1==433 && id2==431)) {
		_c_DsDsS->fill();
	      }
	    }
	  }
	  else if(ncount2==1) {
	    int ipi=0;
	    if(nRes2[111]==1 && nRes2[211]==0 && nRes[-211]==0 )
	      ipi = 111;
	    else if(nRes2[111]==0 && nRes2[211]==1 && nRes[-211]==0 )
	      ipi = 211;
	    else if(nRes2[111]==0 && nRes2[211]==0 && nRes[-211]==1 )
	      ipi =-211;
	    if(ipi==0) continue;
	    matched=true;
	    for(auto const & val : nRes2) {
	      if(val.first==ipi)
		continue;
	      else if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      bool Ddecay = false;
	      Particle mother = p1;
	      while (!mother.parents().empty()) {
		mother = mother.parents()[0];
		if(PID::isCharmMeson(mother.pid()) && mother.pid()!=p1.pid()) {
		  Ddecay = true;
		  break;
		}
	      }
	      mother = p2;
	      while (!mother.parents().empty()) {
		mother = mother.parents()[0];
		if(PID::isCharmMeson(mother.pid()) && mother.pid()!=p1.pid()) {
		  Ddecay = true;
		  break;
		}
	      }
	      if(Ddecay) continue;
	      if((id1==413 || id1==423 ) &&
		 (id2==413 || id2==423 )) {
		_c_DSDSpi->fill();
	      }
	      else if((id1==411 || id1==421 ) &&
		      (id2==413 || id2==423 )) {
		_c_DSDpi->fill();
	      }
	      else if((id1==413 || id1==423 ) &&
		      (id2==411 || id2==421 )) {
		_c_DSDpi->fill();
	      }
	    }
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // R
      Scatter1D R = *_c_hadrons/ *_c_muons;
      double              rval = R.point(0).x();
      pair<double,double> rerr = R.point(0).xErrs();
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_c = _c_DDX->val()*fact;
      double err_c = _c_DDX->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Scatter2D temphisto(refData(6, 1, 2));
      Scatter2DPtr charm;
      book(charm, 6,1,1);
      Scatter2DPtr hadrons;
      book(hadrons, "sigma_hadrons");
      Scatter2DPtr muons;
      book(muons, "sigma_muons"  );
      Scatter2DPtr mult;
      book(mult, 6,1,2);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
	  mult   ->addPoint(x, rval, ex, rerr);
	  hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
	  charm  ->addPoint(x, sig_c, ex, make_pair(err_c,err_c));
	  muons  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	}
	else {
	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	  hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
	  charm  ->addPoint(x, 0., ex, make_pair(0.,.0));
	  muons  ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
      for(unsigned int ix=1;ix<6;++ix) {
	unsigned int imax(0);
	if     (ix<=3) imax = 4;
	else           imax = 3;
	for(unsigned int iy=1;iy<imax;++iy) {
	  double sigma(0),error(0);
	  if(ix==1) {
	    if(iy==1) {
	      sigma = _c_D0D0->val()/picobarn;
	      error = _c_D0D0->err()/picobarn;
	    }
	    else if(iy==2) {
	      sigma = _c_D0D0S->val()/picobarn;
	      error = _c_D0D0S->err()/picobarn;
	    }
	    else if(iy==3) {
	      sigma = _c_D0SD0S->val()/picobarn;
	      error = _c_D0SD0S->err()/picobarn;
	    }
	  }
	  else if(ix==2) {
	    if(iy==1) {
	      sigma = _c_DpDm->val()/picobarn;
	      error = _c_DpDm->err()/picobarn;
	    }
	    else if(iy==2) {
	      sigma = _c_DpDmS->val()/picobarn;
	      error = _c_DpDmS->err()/picobarn;
	    }
	    else if(iy==3) {
	      sigma = _c_DpSDmS->val()/picobarn;
	      error = _c_DpSDmS->err()/picobarn;
	    }
	  }
	  else if(ix==3) {
	    if(iy==1) {
	      sigma = _c_DsDs->val()/picobarn;
	      error = _c_DsDs->err()/picobarn;
	    }
	    else if(iy==2) {
	      sigma = _c_DsDsS->val()/picobarn;
	      error = _c_DsDsS->err()/picobarn;
	    }
	    else if(iy==3) {
	      sigma = _c_DsSDsS->val()/picobarn;
	      error = _c_DsSDsS->err()/picobarn;
	    }
	  }
	  else if(ix==4) {
	    if(iy==1) {
	      sigma = _c_DSDpi->val()/picobarn;
	      error = _c_DSDpi->err()/picobarn;
	    }
	    else if(iy==2) {
	      sigma = _c_DSDSpi->val()/picobarn;
	      error = _c_DSDSpi->err()/picobarn;
	    }
	  }
	  else if(ix==5) {
	    if(iy==1) {
	      sigma = _c_DD->val()/nanobarn;
	      error = _c_DD->err()/nanobarn;
	    }
	    else if(iy==2) {
	      sigma = _c_DDX->val()/nanobarn;
	      error = _c_DDX->err()/nanobarn;
	    }
	  }
	  sigma *= crossSection()/ sumOfWeights();
	  error *= crossSection()/ sumOfWeights();
	  Scatter2D temphisto(refData(ix, 1, iy));
	  Scatter2DPtr mult;
	  book(mult, ix,1,iy);
	  for (size_t b = 0; b < temphisto.numPoints(); b++) {
	    const double x  = temphisto.point(b).x();
	    pair<double,double> ex = temphisto.point(b).xErrs();
	    pair<double,double> ex2 = ex;
	    if(ex2.first ==0.) ex2. first=0.0001;
	    if(ex2.second==0.) ex2.second=0.0001;
	    if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
	      mult   ->addPoint(x, sigma, ex, make_pair(error,error));
	    }
	    else {
	      mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	    }
	  }
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_D0D0, _c_DpDm,_c_DsDs;
    CounterPtr _c_D0D0S, _c_DpDmS,_c_DsDsS;
    CounterPtr _c_D0SD0S, _c_DpSDmS,_c_DsSDsS;
    CounterPtr _c_DD, _c_DDX;
    CounterPtr _c_DSDpi, _c_DSDSpi;
    CounterPtr _c_hadrons, _c_muons;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CLEOC_2008_I777917);


}

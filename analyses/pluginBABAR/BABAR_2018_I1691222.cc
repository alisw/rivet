// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief e+e- > e+e- eta'
  class BABAR_2018_I1691222 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2018_I1691222);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      declare(UnstableParticles(), "UFS");
      // book the histograms
      book(_h_etap,1,1,1);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        } else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    bool findScattered(Particle beam, double& q2) {
      bool found = false;
      Particle scat = beam;
      while (!scat.children().empty()) {
        found = false;
        for (const Particle & p : scat.children()) {
          if (p.pid()==scat.pid()) {
            scat=p;
            found=true;
            break;
          }
        }
        if (!found) break;
      }
      if (!found) return false;
      q2 = -(beam.momentum() - scat.momentum()).mass2();
      return true;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find scattered leptons and calc Q2
      const Beam& beams = apply<Beam>(event, "Beams");
      double q12 = -1, q22 = -1;
      if (!findScattered(beams.beams().first,  q12)) vetoEvent;
      if (!findScattered(beams.beams().second, q22)) vetoEvent;
      if(q22>q12) swap(q12,q22);
      // check the final state
      const FinalState & fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find the meson
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==331)) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched = true;
        for(auto const & val : nRes) {
          if(abs(val.first)==11) {
            if(val.second!=1) {
              matched = false;
              break;
            }
          }
          else if(val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched) {
	  // 2<Q2<10 for both photons bin
	  if(q12>2.&&q12<10.&&q22>2.&&q22<10.) {
	    _h_etap->fill(0,1./sqr(8.));
	  }
	  // 10<Q2<30 for both photons bin
	  else if(q12>10&&q12<30.&&q22>10.&&q22<30.)
	    _h_etap->fill(1,1./sqr(20.));
	  // 10<Q12<30 2<Q22<10
	  else if(q22>2.&&q22<10.&&q12>10.&&q12<30.)
	    _h_etap->fill(2,1./8./20./2.);
	  // 2<Q22<30 30<Q12<60
	  else if(q22>2.&&q22<30.&&q12>30.&&q12<60.)
	    _h_etap->fill(3,1./28./30./2.);
	  // 30<Q2<60 for both photons
	  else if(q22>30.&&q22<60.&&q12>30.&&q12<60.)
	    _h_etap->fill(4,1./sqr(30.));
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_etap, 1e4*crossSection()/femtobarn/sumW());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_etap;
    unsigned int _ncount=0;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2018_I1691222);

}

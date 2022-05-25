// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief BELLE  e+e- > e+e- pi0
  class BELLE_2012_I1114749 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2012_I1114749);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      declare(UnstableParticles(), "UFS");
      // book the histograms
      book(_h_pi0,1,1,1);
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
      double scale = max(q12,q22);
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
      for (const Particle& p : ufs.particles(Cuts::pid==111)) {
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
	  _h_pi0->fill(scale);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize the cross sections
      scale(_h_pi0, crossSection()/femtobarn/sumW());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pi0;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2012_I1114749);

}

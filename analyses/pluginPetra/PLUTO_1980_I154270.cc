// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PLUTO_1980_I154270 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PLUTO_1980_I154270);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      if (fuzzyEquals(sqrtS()/GeV,9.4 ) ||
          fuzzyEquals(sqrtS()/GeV,12.0) ||
          fuzzyEquals(sqrtS()/GeV,13.0) ||
          fuzzyEquals(sqrtS()/GeV,17.0) ||
          fuzzyEquals(sqrtS()/GeV,22.0) ||
          fuzzyEquals(sqrtS()/GeV,27.6) ||
          fuzzyEquals(sqrtS()/GeV,30.2) ||
          fuzzyEquals(sqrtS()/GeV,30.7) ||
          fuzzyEquals(sqrtS()/GeV,31.3)) {
        book(_c_mult, "/TMP/cmult");
        book(_mult, 1, 1, 1);
      }
      else {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      unsigned int nPart(0);
      for (const Particle& p : cfs.particles()) {
        // check if prompt or not
        ConstGenParticlePtr pmother = p.genParticle();
        ConstGenVertexPtr ivertex = pmother->production_vertex();
        bool prompt = true;
        while (ivertex) {
          vector<ConstGenParticlePtr> inparts = HepMCUtils::particles(ivertex, Relatives::PARENTS);
          int n_inparts = inparts.size();
          if (n_inparts < 1) break;
          pmother = inparts[0]; // first mother particle
          int mother_pid = abs(pmother->pdg_id());
          if (mother_pid==PID::K0S || mother_pid==PID::LAMBDA) {
            prompt = false;
            break;
          }
          else if (mother_pid<6) {
            break;
          }
          ivertex = pmother->production_vertex();
        }
        if(prompt) ++nPart;
      }
      _c_mult->fill(sqrtS(), nPart);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = 1./sumOfWeights();
      double val = _c_mult->val()*fact;
      double err = _c_mult->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
        const double x  = temphisto.point(b).x();
        pair<double,double> ex = temphisto.point(b).xErrs();
        pair<double,double> ex2 = ex;
        if(ex2.first ==0.) ex2. first=0.0001;
        if(ex2.second==0.) ex2.second=0.0001;
        if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
          _mult->addPoint(x, val, ex, make_pair(err,err));
        }
        else {
          _mult->addPoint(x, 0., ex, make_pair(0.,.0));
        }
      }
    }

    //@}


  private:

    Profile1DPtr _hist;
    CounterPtr _c_mult;
    Scatter2DPtr _mult;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PLUTO_1980_I154270);


}

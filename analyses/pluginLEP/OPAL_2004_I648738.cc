// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class OPAL_2004_I648738 : public Analysis {
  public:

    /// Constructor
    OPAL_2004_I648738()
      : Analysis("OPAL_2004_I648738"), _sumW(3), _histo_xE(3)
    {    }


    /// @name Analysis methods
    //@{
    void init() {
      declare(FinalState(), "FS");
      declare(ChargedFinalState(), "CFS");
      unsigned int ih=0;
      if (inRange(0.5*sqrtS()/GeV, 4.0, 9.0)) {
	ih = 1;
      }
      else if (inRange(0.5*sqrtS()/GeV, 9.0, 19.0)) {
	ih = 2;
      }
      else if (inRange(0.5*sqrtS()/GeV, 19.0, 30.0)) {
	ih = 3;
      }
      else if (inRange(0.5*sqrtS()/GeV, 45.5, 45.7)) {
	ih = 5;
      }
      else if (inRange(0.5*sqrtS()/GeV, 30.0, 70.0)) {
	ih = 4;
      }
      else if (inRange(0.5*sqrtS()/GeV, 91.5, 104.5)) {
	ih = 6;
      }
      assert(ih>0);
      // book the histograms
      book(_histo_xE[0], ih+5,1,1);
      book(_histo_xE[1], ih+5,1,2);
      if(ih<5) book(_histo_xE[2] ,ih+5,1,3);
      book(_sumW[0], "_sumW_0");
      book(_sumW[1], "_sumW_1");
      book(_sumW[2], "_sumW_2");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the initial quarks/gluons
      Particles initial;
      for (ConstGenParticlePtr p : HepMCUtils::particles(event.genEvent())) {
        ConstGenVertexPtr pv = p->production_vertex();
        const PdgId pid = abs(p->pdg_id());
        if(!( (pid>=1&&pid<=5) || pid ==21) ) continue;
        bool passed = false;
        for (ConstGenParticlePtr pp : HepMCUtils::particles(pv, Relatives::PARENTS)) {
          const PdgId ppid = abs(pp->pdg_id());
          passed = (ppid == PID::ELECTRON || ppid == PID::HIGGS ||
                    ppid == PID::ZBOSON   || ppid == PID::GAMMA);
          if(passed) break;
        }
        if(passed) initial.push_back(Particle(*p));
        }
        if(initial.size()!=2) {
          vetoEvent;
        }
        // type of event
        unsigned int itype=2;
        if(initial[0].pid()==-initial[1].pid()) {
          PdgId pid = abs(initial[0].pid());
          if(pid>=1&&pid<=4)
            itype=0;
          else
            itype=1;
        }
        assert(itype<_histo_xE.size());

      // fill histograms
      _sumW[itype]->fill(2.);
      const Particles& chps = applyProjection<FinalState>(event, "CFS").particles();
      for(const Particle& p : chps) {
        double xE = 2.*p.E()/sqrtS();
	_histo_xE[itype]->fill(xE);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<_histo_xE.size();++ix) {
	if(_sumW[ix]->val()>0.) scale(_histo_xE[ix],1./ *_sumW[ix]);
      }
    }
    //@}


  private:

    vector<CounterPtr> _sumW;

    /// @name Histograms
    //@{
    vector<Histo1DPtr> _histo_xE;
    //@}


  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2004_I648738);

}

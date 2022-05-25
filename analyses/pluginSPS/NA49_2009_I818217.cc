// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  class NA49_2009_I818217 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{
    /// Constructor
    NA49_2009_I818217() : Analysis("NA49_2009_I818217"){}

    void init() {
      declare(FinalState(), "FS");
      book(_h_dndxf, 1,1,1);
      book(_p_mptxf, 2,1,1);
      book(_h_dndy,  3,1,1);
      book(_c_ninel, "nInelastic");
    }

    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      const float SRT = event.sqrtS();

      // Inelastic events selection
      if (numParticles <= 2) {
        MSG_DEBUG("Elastic event");
        vetoEvent;
      }
      _c_ninel -> fill();

      // Plot distributions
      for(const Particle& p : fs.particles()) {
        if(p.pid() == PID::PROTON){
          double xF = p.pz() / (SRT/2.);
          _h_dndxf -> fill(xF);
          _h_dndy  -> fill(p.rapidity());
          _p_mptxf -> fill(xF, p.pt());
        }
      }
    }

    void finalize() {
      scale(_h_dndxf, 1./ *_c_ninel);   // Scale by the number of inelastic events
      vector<YODA::HistoBin1D>& bins = _h_dndxf -> bins(); // Get histogram bins
      for (auto b : bins) b.scaleW(1./b.xWidth());  // Scale by the bin width (dxF)

      scale(_h_dndy, 1./ *_c_ninel);
      vector<YODA::HistoBin1D>& binsy = _h_dndy -> bins();
      for (auto by : binsy) by.scaleW(1./by.xWidth());
    }

  private:
    CounterPtr   _c_ninel;  // Counter of inelastic events
    Histo1DPtr   _h_dndxf;  // dn/dxf histogram
    Histo1DPtr   _h_dndy;   // dn/dy histogram
    Profile1DPtr _p_mptxf;  // mean pT vs xF profile
  };

  RIVET_DECLARE_PLUGIN(NA49_2009_I818217);
}

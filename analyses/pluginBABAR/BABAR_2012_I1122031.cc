// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Xs gamma
  class BABAR_2012_I1122031 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1122031);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h, 1, 1, 1);
      book(_nBottom, "TMP/BottomCounter");
    }

    void findDecayProducts(const Particle& mother,
                           unsigned int& nK0, unsigned int& nKp, unsigned int& nKm,
                           FourMomentum& ptot) {
      for (const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS ) {
          ++nKp;
          ptot += p.momentum();
        }
        else if (id == PID::KMINUS ) {
          ++nKm;
          ptot += p.momentum();
        }
        else if (id == PID::K0S) {
          ++nK0;
          ptot += p.momentum();
        }
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
          ptot += p.momentum();
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nK0, nKp, nKm, ptot);
        }
        else
          ptot += p.momentum();
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over bottoms
      for (const Particle& bottom : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==521 ||
										     Cuts::abspid==511)) {
	// remove mixing entries etc
	for (const Particle & child : bottom.children())
          if (child.abspid() == 511 || child.pid()==bottom.pid() ) continue;
        _nBottom->fill();
	FourMomentum pgamma(0.,0.,0.,0.);
	unsigned int ngamma = 0;
        for (const Particle & child : bottom.children()) {
	  if (child.pid() == PID::PHOTON) {
            ngamma += 1;
            pgamma += child.momentum();
          }
	}
	if (ngamma != 1) continue;
        unsigned int nK0(0),nKp(0),nKm(0);
        FourMomentum p_tot(0,0,0,0);
        findDecayProducts(bottom, nK0, nKp, nKm, p_tot);
        unsigned int nk = nKp-nKm+nK0;
        if (nk % 2 == 1) {
          p_tot -= pgamma;
          _h->fill(p_tot.mass()/GeV);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1e6/_nBottom->sumW());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _nBottom;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1122031);

}

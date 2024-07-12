// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> s gamma
  class BELLE_2022_I2167323 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2167323);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==521 || Cuts::abspid==511), "UFS");
      // Book histograms
      book(_h_spectrum, 1, 1, 1);
      book(_nBottom, "TMP/BottomCounter");
    }

    void findDecayProducts(const Particle& mother,
                           unsigned int& nK0, unsigned int& nKp,
			   unsigned int& nKm) {
      for (const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS )      ++nKp;
        else if (id == PID::KMINUS ) ++nKm;
        else if (id == PID::K0S)     ++nK0;
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
	  continue;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nK0, nKp, nKm);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over bottoms
      for (const Particle& bottom : apply<UnstableParticles>(event, "UFS").particles()) {
	// remove mixing entries etc
	if(bottom.children()[0].abspid()==bottom.abspid()) continue;
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
        findDecayProducts(bottom, nK0, nKp, nKm);
        unsigned int nk = nKp-nKm+nK0;
        if (nk % 2 == 1) {
	  const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(bottom.momentum().betaVec());
	  double eGamma = boost.transform(pgamma).E();
	  _h_spectrum->fill(eGamma);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spectrum, 1e4/_nBottom->sumW());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spectrum;
    CounterPtr _nBottom;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2167323);

}

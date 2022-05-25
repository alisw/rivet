// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALEPH_2014_I1267648 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_2014_I1267648);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_pip0 , 1, 1, 1);
      book(_h_pi2p0, 2, 1, 1);
      book(_h_pi3p0, 3, 1, 1);
      book(_h_3pi  , 4, 1, 1);
      book(_h_3pip0, 5, 1, 1);
    }


    void findDecayProducts(const Particle &mother, unsigned int &nstable, unsigned int &npip,
                           unsigned int &npim, unsigned int &npi0, FourMomentum &ptot) {
      for (const Particle &p : mother.children()) {
        int id = p.pid();
        if (id == PID::KPLUS || id == PID::KMINUS) {
          ++nstable;
          ptot += p.momentum();
        }
        else if (id == PID::PIPLUS) {
          ++npip;
          ++nstable;
          ptot += p.momentum();
        }
        else if (id == PID::PIMINUS) {
          ++npim;
          ++nstable;
          ptot += p.momentum();
        }
        else if (id == PID::PI0) {
          ++nstable;
          ++npi0;
          ptot += p.momentum();
        }
        else if (id == PID::PHOTON)  continue;
        else if (!p.children().empty())  findDecayProducts(p, nstable, npip, npim, npi0, ptot);
        else  ++nstable;
      }
    }
    

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over taus
      for (const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
        FourMomentum ptot;
        unsigned int nstable(0), npip(0), npim(0), npi0(0);
        findDecayProducts(tau,nstable,npip,npim,npi0,ptot);
        // tau -> pi pi0 nu_tau (both charges)
        if (npim==1 && npi0==1 && nstable==3)  _h_pip0->fill(ptot.mass2());
        // tau -> pi pi0 pi0 nu_tau (both charges)
        else if (npim==1 && npi0==2 && nstable==4)  _h_pi2p0->fill(ptot.mass2());
        //    tau -> pi pi0 pi0 pi0         (3,1,1)
        else if (npim==1 && npi0==3 && nstable==5)  _h_pi3p0->fill(ptot.mass2());
        //    tau -> 3 charged pions        (4,1,1)
        else if (npim==2 && npip==1 && nstable==4)  _h_3pi->fill(ptot.mass2());
        //    tau -> 3 charged pions + pi0  (5,1,1)
        else if (npim==2 && npip==1 && npi0==1 && nstable==5)  _h_3pip0->fill(ptot.mass2());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_pip0);  // normalize to unity
      normalize(_h_pi2p0); // normalize to unity
      normalize(_h_pi3p0); // normalize to unity
      normalize(_h_3pi);   // normalize to unity
      normalize(_h_3pip0); // normalize to unity

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_pip0;
    Histo1DPtr _h_pi2p0;
    Histo1DPtr _h_pi3p0;
    Histo1DPtr _h_3pi;
    Histo1DPtr _h_3pip0;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALEPH_2014_I1267648);

}

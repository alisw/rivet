// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Production of light-flavor hadrons in pp collisions at 13 TeV
  class ALICE_2020_I1797443 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2020_I1797443);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // All final-state and unstable particles within
      // the given rapidity acceptance
      const UnstableParticles up(Cuts::absrap < 0.5);
      declare(up, "up");

      // Book histograms
      book(_h["pi"], 1, 1, 1);
      book(_h["k"], 2, 1, 1);
      book(_h["k0s"], 3, 1, 1);
      book(_h["k*0"], 4, 1, 1);
      book(_h["phi"], 5, 1, 1);
      book(_h["p"], 6, 1, 1);
      book(_h["l0"], 7, 1, 1);
      book(_h["xi"], 8, 1, 1);
      book(_h["omega"], 9, 1, 1);
      // Book ratios
      book(_s["k/pi"], 21, 1, 1);
      book(_s["k*0/pi"], 22, 1, 1);
      book(_s["k0s/pi"], 23, 1, 1);
      book(_s["phi/pi"], 24, 1, 1);
      book(_s["p/pi"], 42, 1, 1);
      book(_s["l0/k0s"], 43, 1, 1);
      book(_s["omega/phi"], 44, 1, 1);
      book(_s["xi/phi"], 45, 1, 1);
      // Book temporary histograms for ratios
      book(_h["k_for_k/pi"],  "TMP/k_for_k/pi",  refData(21, 1, 1));
      book(_h["pi_for_k/pi"], "TMP/pi_for_k/pi", refData(21, 1, 1));
      book(_h["k*0_for_k*0/pi"], "TMP/k*0_for_k*0/pi", refData(22, 1, 1));
      book(_h["pi_for_k*0/pi"],  "TMP/pi_for_k*0/pi",  refData(22, 1, 1));
      book(_h["k0s_for_k0s/pi"], "TMP/k0s_for_k0s/pi", refData(23, 1, 1));
      book(_h["pi_for_k0s/pi"],  "TMP/pi_for_k0s/pi",  refData(23, 1, 1));
      book(_h["phi_for_phi/pi"], "TMP/phi_for_phi/pi", refData(24, 1, 1));
      book(_h["pi_for_phi/pi"],  "TMP/pi_for_phi/pi",  refData(24, 1, 1));
      book(_h["p_for_p/pi"],  "TMP/p_for_p/pi",  refData(42, 1, 1));
      book(_h["pi_for_p/pi"], "TMP/pi_for_p/pi", refData(42, 1, 1));
      book(_h["l0_for_l0/k0s"],  "TMP/l0_for_l0/k0s",  refData(43, 1, 1));
      book(_h["k0s_for_l0/k0s"], "TMP/k0s_for_l0/k0s", refData(43, 1, 1));
      book(_h["omega_for_omega/phi"], "TMP/omega_for_omega/phi", refData(44, 1, 1));
      book(_h["phi_for_omega/phi"],   "TMP/phi_for_omega/phi",   refData(44, 1, 1));
      book(_h["xi_for_xi/phi"],  "TMP/xi_for_xi/phi",  refData(45, 1, 1));
      book(_h["phi_for_xi/phi"], "TMP/phi_for_xi/phi", refData(45, 1, 1));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& up = apply<UnstableParticles>(event, "up").particles();

      for (const Particle& p : up) {
        if (p.abspid() == 211) {
          _h["pi"]->fill(p.pT() / GeV);
          _h["pi_for_k/pi"]->fill(p.pT() / GeV);
          _h["pi_for_k*0/pi"]->fill(p.pT() / GeV);
          _h["pi_for_k0s/pi"]->fill(p.pT() / GeV);
          _h["pi_for_phi/pi"]->fill(p.pT() / GeV);
          _h["pi_for_p/pi"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 321) {
          _h["k"]->fill(p.pT() / GeV);
          _h["k_for_k/pi"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 310) {
          _h["k0s"]->fill(p.pT() / GeV);
          _h["k0s_for_k0s/pi"]->fill(p.pT() / GeV);
          _h["k0s_for_l0/k0s"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 2212) {
          _h["p"]->fill(p.pT() / GeV);
          _h["p_for_p/pi"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 313) {
          _h["k*0"]->fill(p.pT() / GeV);
          _h["k*0_for_k*0/pi"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 333) {
          _h["phi"]->fill(p.pT() / GeV);
          _h["phi_for_phi/pi"]->fill(p.pT() / GeV);
          _h["phi_for_omega/phi"]->fill(p.pT() / GeV);
          _h["phi_for_xi/phi"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 3122) {
          _h["l0"]->fill(p.pT() / GeV);
          _h["l0_for_l0/k0s"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 3312) {
          _h["xi"]->fill(p.pT() / GeV);
          _h["xi_for_xi/phi"]->fill(p.pT() / GeV);
        }
        else if (p.abspid() == 3334) {
          _h["omega"]->fill(p.pT() / GeV);
          _h["omega_for_omega/phi"]->fill(p.pT() / GeV);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h, 1./sumOfWeights());
      scale({_h["k0s_for_k0s/pi"], _h["phi_for_phi/pi"], _h["k0s_for_l0/k0s"]}, 2.);

      divide(_h["k_for_k/pi"], _h["pi_for_k/pi"], _s["k/pi"]);
      divide(_h["k*0_for_k*0/pi"], _h["pi_for_k*0/pi"], _s["k*0/pi"]);
      divide(_h["k0s_for_k0s/pi"], _h["pi_for_k0s/pi"], _s["k0s/pi"]);
      divide(_h["phi_for_phi/pi"], _h["pi_for_phi/pi"], _s["phi/pi"]);
      divide(_h["p_for_p/pi"], _h["pi_for_p/pi"], _s["p/pi"]);
      divide(_h["l0_for_l0/k0s"], _h["k0s_for_l0/k0s"], _s["l0/k0s"]);
      divide(_h["omega_for_omega/phi"], _h["phi_for_omega/phi"], _s["omega/phi"]);
      divide(_h["xi_for_xi/phi"], _h["phi_for_xi/phi"], _s["xi/phi"]);

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr> _s;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2020_I1797443);

}

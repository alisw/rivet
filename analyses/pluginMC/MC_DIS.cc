// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

class MC_DIS : public Analysis {

public:

  /// Constructor
  DEFAULT_RIVET_ANALYSIS_CTOR(MC_DIS);

  /// @name Analysis methods
  /// @{

  /// Book histograms and initialise projections before the run
  void init() {

    // Initialise and register projections. Note that the definition
    // of the scattered lepton can be influenced by specifying
    // options as declared in the .info file.
    declare(FinalState(),  "FS");
    DISLepton lepton(options());
    declare(lepton, "Lepton");
    declare(DISKinematics(lepton), "Kinematics");

    book(_h_charge_electron,    "chargeelectron", 2, -1.0, 1.0);

    vector<double> bin_edges_of_x;
    for (size_t i = 0; i < 100 + 1; i++) bin_edges_of_x.push_back(0.000001*pow(1.0/0.000001, i/100.0));

    book(_h_x, "x",  bin_edges_of_x);
    book(_h_eminuspz, "eminuspz", 240, 0.0, 60.0);
    book(_h_etot_remnant, "etotremnant", 100, 0.0, 1000.0);
    book(_h_pt_remnant, "ptremnant", 50, 0.0, 5.0);

    book(_h_pttot, "pttot", 200, 0.0, 200.0);
    book(_h_pttot_leptons, "pttotleptons", 200, 0.0, 200.0);
    book(_h_pttot_hadrons, "pttothadrons", 200, 0.0, 200.0);
    book(_h_pttot_gamma, "pttotgamma", 200, 0.0, 200.0);

    book(_h_e_electron, "eelectron",240, 0.0, 60.0);
    book(_h_pt_electron, "ptelectron", 240, 0.0, 60.0);
    book(_h_y, "y", 100, 0.0, 1.0);
    book(_h_W2, "W2", 100, 0.0, 100000);

    vector<double> bin_edges_of_Q2;
    for (size_t i = 0; i < 100 + 1; i++) bin_edges_of_Q2.push_back(0.1*pow(100000.0/0.1, i/100.0));

    book(_h_Q2, "Q2", bin_edges_of_Q2);
    book(_h_gammahad, "gammahad", 180, 0.0, 180);
    book(_h_theta_electron, "thetaelectron", 180, 0.0, 180);
  }

  /// Perform the per-event analysis
  void analyze(const Event& event) {

    /// We analyze event an extract DIS kinematics
    const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
    const DISLepton& dl = apply<DISLepton>(event,"Lepton");

    const double q2 = dk.Q2();
    const double x = dk.x();
    const double y = dk.y();
    const double W2 = dk.W2();
    const double gammahad = dk.gammahad()/degree;

    _h_x->fill(x);
    _h_y->fill(y);
    _h_W2->fill(W2);
    _h_Q2->fill(q2);
    _h_gammahad->fill(gammahad);
    _h_theta_electron->fill(dl.out().angle(dk.beamHadron().mom())/degree);
    _h_e_electron->fill(dl.out().E());
    _h_pt_electron->fill(dl.out().pT());
    _h_charge_electron->fill(0.5*(dl.in().charge() > 0 ? 1. : -1));

    double eminuspz = 0;
    double etot_remnant = 0;
    double pttot = 0; /// transverse momentum of all particles but the scattered lepton
    double pttot_leptons = 0; /// transverse momentum of all leptons but the scattered one
    double pttot_hadrons = 0; /// transverse momentum of all hadrons
    double pttot_gamma = 0;   /// transverse momentum of all gammas
    const FinalState& fs = apply<FinalState>(event, "FS");
    for (const Particle& p: fs.particles()) {
      eminuspz += ( p.E() + p.pz()*dl.pzSign());
      if ( p.genParticle() == dl.out().genParticle() ) continue;
      pttot += p.pT();
      if ( p.isLepton() ) pttot_leptons += p.pT();
      if ( p.abspid() == PID::PHOTON ) pttot_gamma += p.pT();
      if ( p.isVisible() && !p.isLepton() && !(p.abspid() == PID::PHOTON) ) pttot_hadrons += p.pT();

      if ( p.abseta() < 6 ) continue;
      etot_remnant += p.E();
      _h_pt_remnant->fill(p.pT());
    }
    _h_eminuspz->fill(eminuspz);
    _h_etot_remnant->fill(etot_remnant);
    _h_pttot->fill(pttot);
    _h_pttot_leptons->fill(pttot_leptons);
    _h_pttot_hadrons->fill(pttot_hadrons);
    _h_pttot_gamma->fill(pttot_gamma);
  }

  /// Normalise histograms etc., after the run
  void finalize() {
    scale(_h_charge_electron, crossSection()/sumOfWeights());
    normalize({_h_y, _h_W2, _h_x, _h_Q2, _h_gammahad,
               _h_eminuspz,
               _h_pt_remnant,
               _h_etot_remnant,
               _h_pttot, _h_pttot_leptons, _h_pttot_hadrons, _h_pttot_gamma,
               _h_e_electron, _h_pt_electron, _h_theta_electron});

  }
  /// @}

private:

  /// @name Histograms
  /// @{
  Histo1DPtr _h_charge_electron;
  Histo1DPtr _h_y, _h_W2, _h_x, _h_Q2, _h_gammahad;
  Histo1DPtr _h_eminuspz;
  Histo1DPtr _h_pt_remnant;
  Histo1DPtr _h_etot_remnant;
  Histo1DPtr _h_pttot, _h_pttot_leptons, _h_pttot_hadrons, _h_pttot_gamma;
  Histo1DPtr _h_e_electron, _h_pt_electron, _h_theta_electron;

  /// @}
};



DECLARE_RIVET_PLUGIN(MC_DIS);

}

// -*- C++ -*-
#ifndef RIVET_SmearingFunctions_HH
#define RIVET_SmearingFunctions_HH

#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"
#include <random>

namespace Rivet {


  /// Return a uniformly sampled random number between 0 and 1
  /// @todo Move to (math?)utils
  /// @todo Need to isolate random generators to a single thread
  inline double rand01() {
    //return rand() / (double)RAND_MAX;
    static random_device rd;
    static mt19937 gen(rd());
    return generate_canonical<double, 10>(gen);
  }


  /// @name General particle & momentum efficiency and smearing functions
  //@{

  /// Take a Particle and return 0
  inline double PARTICLE_FN0(const Particle& p) { return 0; }
  /// Take a Particle and return 1
  inline double PARTICLE_FN1(const Particle& p) { return 1; }
  /// Take a Particle and return it unmodified
  inline Particle PARTICLE_SMEAR_IDENTITY(const Particle& p) { return p; }


  /// Take a FourMomentum and return 0
  inline double P4_FN0(const FourMomentum& p) { return 0; }
  /// Take a FourMomentum and return 1
  inline double P4_FN1(const FourMomentum& p) { return 1; }
  /// Take a FourMomentum and return it unmodified
  inline FourMomentum P4_SMEAR_IDENTITY(const FourMomentum& p) { return p; }

  /// Smear a FourMomentum's energy using a Gaussian of absolute width @a resolution
  /// @todo Also make jet versions that update/smear constituents?
  inline FourMomentum P4_SMEAR_E_GAUSS(const FourMomentum& p, double resolution) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());
    normal_distribution<> d(p.E(), resolution);
    const double mass = p.mass2() > 0 ? p.mass() : 0; //< numerical carefulness...
    const double smeared_E = max(d(gen), mass); //< can't let the energy go below the mass!
    return FourMomentum::mkEtaPhiME(p.eta(), p.phi(), mass, smeared_E);
  }

  /// Smear a FourMomentum's transverse momentum using a Gaussian of absolute width @a resolution
  inline FourMomentum P4_SMEAR_PT_GAUSS(const FourMomentum& p, double resolution) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());
    normal_distribution<> d(p.pT(), resolution);
    const double smeared_pt = max(d(gen), 0.);
    const double mass = p.mass2() > 0 ? p.mass() : 0; //< numerical carefulness...
    return FourMomentum::mkEtaPhiMPt(p.eta(), p.phi(), mass, smeared_pt);
  }

  /// Smear a FourMomentum's mass using a Gaussian of absolute width @a resolution
  inline FourMomentum P4_SMEAR_MASS_GAUSS(const FourMomentum& p, double resolution) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());
    normal_distribution<> d(p.mass(), resolution);
    const double smeared_mass = max(d(gen), 0.);
    return FourMomentum::mkEtaPhiMPt(p.eta(), p.phi(), smeared_mass, p.pT());
  }


  /// Take a Vector3 and return 0
  inline double P3_FN0(const Vector3& p) { return 0; }
  /// Take a Vector3 and return 1
  inline double P3_FN1(const Vector3& p) { return 1; }
  /// Take a Vector3 and return it unmodified
  inline Vector3 P3_SMEAR_IDENTITY(const Vector3& p) { return p; }

  /// Smear a Vector3's length using a Gaussian of absolute width @a resolution
  inline Vector3 P3_SMEAR_LEN_GAUSS(const Vector3& p, double resolution) {
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());
    normal_distribution<> d(p.mod(), resolution);
    const double smeared_mod = max(d(gen), 0.); //< can't let the energy go below the mass!
    return smeared_mod * p.unit();
  }

  //@}


  /// @name Electron efficiency and smearing functions
  //@{

  /// ATLAS Run 1 electron tracking efficiency
  /// @todo How to use this in combination with reco eff?
  inline double ELECTRON_TRKEFF_ATLAS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 0.1*GeV) return 0;
    if (e.abseta() < 1.5) {
      if (e.pT() < 1*GeV) return 0.73;
      if (e.pT() < 100*GeV) return 0.95;
      return 0.99;
    } else {
      if (e.pT() < 1*GeV) return 0.50;
      if (e.pT() < 100*GeV) return 0.83;
      else return 0.90;
    }
  }


  /// ATLAS Run 2 electron tracking efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_TRKEFF_ATLAS_RUN2(const Particle& e) {
    return ELECTRON_TRKEFF_ATLAS_RUN1(e);
  }

  /// ATLAS Run 1 electron reconstruction efficiency
  /// @todo Include reco eff (but no e/y discrimination) in forward region
  /// @todo How to use this in combination with tracking eff?
  inline double ELECTRON_EFF_ATLAS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 10*GeV) return 0;
    return (e.abseta() < 1.5) ? 0.95 : 0.85;
  }


  /// ATLAS Run 2 electron reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_EFF_ATLAS_RUN2(const Particle& e) {
    return ELECTRON_EFF_ATLAS_RUN1(e);
  }


  /// ATLAS Run 1 electron reco smearing
  inline Particle ELECTRON_SMEAR_ATLAS_RUN1(const Particle& e) {
    static const vector<double> edges_eta = {0., 2.5, 3.};
    static const vector<double> edges_pt = {0., 0.1, 25.};
    static const vector<double> e2s = {0.000, 0.015, 0.005,
                                       0.005, 0.005, 0.005,
                                       0.107, 0.107, 0.107};
    static const vector<double> es = {0.00, 0.00, 0.05,
                                      0.05, 0.05, 0.05,
                                      2.08, 2.08, 2.08};
    static const vector<double> cs = {0.00, 0.00, 0.25,
                                      0.25, 0.25, 0.25,
                                      0.00, 0.00, 0.00};

    const int i_eta = binIndex(e.abseta(), edges_eta, true);
    const int i_pt = binIndex(e.pT()/GeV, edges_pt, true);
    const int i = i_eta*edges_pt.size() + i_pt;

    // Calculate absolute resolution in GeV
    const double c1 = sqr(e2s[i]), c2 = sqr(es[i]), c3 = sqr(cs[i]);
    const double resolution = sqrt(c1*e.E2() + c2*e.E() + c3) * GeV;

    // normal_distribution<> d(e.E(), resolution);
    // const double mass = e.mass2() > 0 ? e.mass() : 0; //< numerical carefulness...
    // const double smeared_E = max(d(gen), mass); //< can't let the energy go below the mass!
    // return Particle(e.pid(), FourMomentum::mkEtaPhiME(e.eta(), e.phi(), mass, smeared_E));
    return Particle(e.pid(), P4_SMEAR_E_GAUSS(e, resolution));
  }


  /// ATLAS Run 2 electron reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle ELECTRON_SMEAR_ATLAS_RUN2(const Particle& e) {
    return ELECTRON_SMEAR_ATLAS_RUN1(e);
  }


  /// @todo Add charge flip efficiency?


  /// CMS Run 1 electron tracking efficiency
  /// @todo How to use this in combination with reco eff?
  inline double ELECTRON_TRKEFF_CMS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 0.1*GeV) return 0;
    if (e.abseta() < 1.5) {
      return (e.pT() < 1*GeV) ? 0.70 : 0.95;
    } else {
      return (e.pT() < 1*GeV) ? 0.60 : 0.85;
    }
  }


  /// CMS Run 2 electron tracking efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_TRKEFF_CMS_RUN2(const Particle& e) {
    return ELECTRON_TRKEFF_CMS_RUN1(e);
  }


  /// CMS Run 1 electron reconstruction efficiency
  /// @todo How to use this in combination with tracking eff?
  inline double ELECTRON_EFF_CMS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 10*GeV) return 0;
    return (e.abseta() < 1.5) ? 0.95 : 0.85;
  }


  /// CMS Run 2 electron reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_EFF_CMS_RUN2(const Particle& e) {
    return ELECTRON_EFF_CMS_RUN1(e);
  }


  /// @brief CMS electron energy smearing, preserving direction
  ///
  /// Calculate resolution
  /// for pT > 0.1 GeV, E resolution = |eta| < 0.5 -> sqrt(0.06^2 + pt^2 * 1.3e-3^2)
  ///                                  |eta| < 1.5 -> sqrt(0.10^2 + pt^2 * 1.7e-3^2)
  ///                                  |eta| < 2.5 -> sqrt(0.25^2 + pt^2 * 3.1e-3^2)
  inline Particle ELECTRON_SMEAR_CMS_RUN1(const Particle& e) {
    // Calculate absolute resolution in GeV from functional form
    double resolution = 0;
    const double abseta = e.abseta();
    if (e.pT() > 0.1*GeV && abseta < 2.5) { //< should be a given from efficiencies
      if (abseta < 0.5) {
        resolution = add_quad(0.06, 1.3e-3 * e.pT()/GeV) * GeV;
      } else if (abseta < 1.5) {
        resolution = add_quad(0.10, 1.7e-3 * e.pT()/GeV) * GeV;
      } else { // still |eta| < 2.5
        resolution = add_quad(0.25, 3.1e-3 * e.pT()/GeV) * GeV;
      }
    }

    // normal_distribution<> d(e.E(), resolution);
    // const double mass = e.mass2() > 0 ? e.mass() : 0; //< numerical carefulness...
    // const double smeared_E = max(d(gen), mass); //< can't let the energy go below the mass!
    // return Particle(e.pid(), FourMomentum::mkEtaPhiME(e.eta(), e.phi(), mass, smeared_E));
    return Particle(e.pid(), P4_SMEAR_E_GAUSS(e, resolution));
  }


  /// CMS Run 2 electron reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle ELECTRON_SMEAR_CMS_RUN2(const Particle& e) {
    return ELECTRON_SMEAR_CMS_RUN1(e);
  }

  //@}



  /// @name Photon efficiency and smearing functions
  //@{

  /// @todo Photon efficiency and smearing

  //@}



  /// @name Muon efficiency and smearing functions
  //@{

  /// ATLAS Run 1 muon tracking efficiency
  /// @todo How to use this in combination with reco eff?
  inline double MUON_TRKEFF_ATLAS_RUN1(const Particle& m) {
    if (m.abseta() > 2.5) return 0;
    if (m.pT() < 0.1*GeV) return 0;
    if (m.abseta() < 1.5) {
      return (m.pT() < 1*GeV) ? 0.75 : 0.99;
    } else {
      return (m.pT() < 1*GeV) ? 0.70 : 0.98;
    }
  }

  /// ATLAS Run 2 muon tracking efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double MUON_TRKEFF_ATLAS_RUN2(const Particle& m) {
    return MUON_TRKEFF_ATLAS_RUN1(m);
  }


  /// ATLAS Run 1 muon reco efficiency
  /// @todo How to use this in combination with tracking eff?
  inline double MUON_EFF_ATLAS_RUN1(const Particle& m) {
    if (m.abseta() > 2.7) return 0;
    if (m.pT() < 10*GeV) return 0;
    return (m.abseta() < 1.5) ? 0.95 : 0.85;
  }

  /// ATLAS Run 2 muon reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double MUON_EFF_ATLAS_RUN2(const Particle& m) {
    return MUON_EFF_ATLAS_RUN1(m);
  }


  /// ATLAS Run 1 muon reco smearing
  inline Particle MUON_SMEAR_ATLAS_RUN1(const Particle& m) {
    static const vector<double> edges_eta = {0, 1.5, 2.5};
    static const vector<double> edges_pt = {0, 0.1, 1.0, 10., 200.};
    static const vector<double> res = {0., 0.03, 0.02, 0.03, 0.05,
                                       0., 0.04, 0.03, 0.04, 0.05};

    const int i_eta = binIndex(m.abseta(), edges_eta, true);
    const int i_pt = binIndex(m.pT()/GeV, edges_pt, true);
    const int i = i_eta*edges_pt.size() + i_pt;

    const double resolution = res[i];

    // Smear by a Gaussian centered on the current pT, with width given by the resolution
    // normal_distribution<> d(m.pT(), resolution*m.pT());
    // const double smeared_pt = max(d(gen), 0.);
    // const double mass = m.mass2() > 0 ? m.mass() : 0; //< numerical carefulness...
    // return Particle(m.pid(), FourMomentum::mkEtaPhiMPt(m.eta(), m.phi(), mass, smeared_pt));
    return Particle(m.pid(), P4_SMEAR_PT_GAUSS(m, resolution*m.pT()));
  }

  /// ATLAS Run 2 muon reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle MUON_SMEAR_ATLAS_RUN2(const Particle& m) {
    return MUON_SMEAR_ATLAS_RUN1(m);
  }


  /// CMS Run 1 muon tracking efficiency
  /// @todo How to use this in combination with reco eff?
  /// @note Eff values currently identical to those in ATLAS (AB, 2016-04-12)
  inline double MUON_TRKEFF_CMS_RUN1(const Particle& m) {
    if (m.abseta() > 2.5) return 0;
    if (m.pT() < 0.1*GeV) return 0;
    if (m.abseta() < 1.5) {
      return (m.pT() < 1*GeV) ? 0.75 : 0.99;
    } else {
      return (m.pT() < 1*GeV) ? 0.70 : 0.98;
    }
  }

  /// CMS Run 2 muon tracking efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double MUON_TRKEFF_CMS_RUN2(const Particle& m) {
    return MUON_TRKEFF_CMS_RUN1(m);
  }


  /// CMS Run 1 muon reco efficiency
  /// @todo How to use this in combination with tracking eff?
  inline double MUON_EFF_CMS_RUN1(const Particle& m) {
    if (m.abseta() > 2.4) return 0;
    if (m.pT() < 10*GeV) return 0;
    return 0.95 * (m.abseta() < 1.5 ? 1 : exp(0.5 - 5e-4*m.pT()/GeV));
  }

  /// CMS Run 2 muon reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double MUON_EFF_CMS_RUN2(const Particle& m) {
    return MUON_EFF_CMS_RUN1(m);
  }


  /// CMS Run 1 muon reco smearing
  inline Particle MUON_SMEAR_CMS_RUN1(const Particle& m) {
    // Calculate fractional resolution
    // for pT > 0.1 GeV, mom resolution = |eta| < 0.5 -> sqrt(0.01^2 + pt^2 * 2.0e-4^2)
    //                                    |eta| < 1.5 -> sqrt(0.02^2 + pt^2 * 3.0e-4^2)
    //                                    |eta| < 2.5 -> sqrt(0.05^2 + pt^2 * 2.6e-4^2)
    double resolution = 0;
    const double abseta = m.abseta();
    if (m.pT() > 0.1*GeV && abseta < 2.5) {
      if (abseta < 0.5) {
        resolution = add_quad(0.01, 2.0e-4 * m.pT()/GeV);
      } else if (abseta < 1.5) {
        resolution = add_quad(0.02, 3.0e-4 * m.pT()/GeV);
      } else { // still |eta| < 2.5... but isn't CMS' mu acceptance < 2.4?
        resolution = add_quad(0.05, 2.6e-4 * m.pT()/GeV);
      }
    }

    // Smear by a Gaussian centered on the current pT, with width given by the resolution
    // normal_distribution<> d(m.pT(), resolution*m.pT());
    // const double smeared_pt = max(d(gen), 0.);
    // const double mass = m.mass2() > 0 ? m.mass() : 0; //< numerical carefulness...
    // return Particle(m.pid(), FourMomentum::mkEtaPhiMPt(m.eta(), m.phi(), mass, smeared_pt));
    return Particle(m.pid(), P4_SMEAR_PT_GAUSS(m, resolution*m.pT()));
  }

  /// CMS Run 2 muon reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle MUON_SMEAR_CMS_RUN2(const Particle& m) {
    return MUON_SMEAR_CMS_RUN1(m);
  }

  //@}



  /// @name Tau efficiency and smearing functions
  //@{

  /// @brief ATLAS Run 1 8 TeV tau efficiencies (medium working point)
  ///
  /// Taken from http://arxiv.org/pdf/1412.7086.pdf
  ///   20-40 GeV 1-prong LMT eff|mis = 0.66|1/10, 0.56|1/20, 0.36|1/80
  ///   20-40 GeV 3-prong LMT eff|mis = 0.45|1/60, 0.38|1/100, 0.27|1/300
  ///   > 40 GeV 1-prong LMT eff|mis = 0.66|1/15, 0.56|1/25, 0.36|1/80
  ///   > 40 GeV 3-prong LMT eff|mis = 0.45|1/250, 0.38|1/400, 0.27|1/1300
  inline double TAU_EFF_ATLAS_RUN1(const Particle& t) {
    if (t.abseta() > 2.5) return 0; //< hmm... mostly
    double pThadvis = 0;
    Particles chargedhadrons;
    for (const Particle& p : t.children()) {
      if (p.isHadron()) {
        pThadvis += p.pT(); //< right definition? Paper is unclear
        if (p.charge3() != 0 && p.abseta() < 2.5 && p.pT() > 1*GeV) chargedhadrons += p;
      }
    }
    if (chargedhadrons.empty()) return 0; //< leptonic tau
    if (pThadvis < 20*GeV) return 0; //< below threshold
    if (pThadvis < 40*GeV) {
      if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.56 : 1/20.;
      if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.38 : 1/100.;
    } else {
      if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.56 : 1/25.;
      if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.38 : 1/400.;
    }
    return 0;
  }


  /// @brief ATLAS Run 2 13 TeV tau efficiencies (medium working point)
  ///
  /// From https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PUBNOTES/ATL-PHYS-PUB-2015-045/ATL-PHYS-PUB-2015-045.pdf
  ///   LMT 1 prong efficiency/mistag = 0.6|1/30, 0.55|1/50, 0.45|1/120
  ///   LMT 3 prong efficiency/mistag = 0.5|1/30, 0.4|1/110, 0.3|1/300
  inline double TAU_EFF_ATLAS_RUN2(const Particle& t) {
    if (t.abseta() > 2.5) return 0; //< hmm... mostly
    double pThadvis = 0;
    Particles chargedhadrons;
    for (const Particle& p : t.children()) {
      if (p.isHadron()) {
        pThadvis += p.pT(); //< right definition? Paper is unclear
        if (p.charge3() != 0 && p.abseta() < 2.5 && p.pT() > 1*GeV) chargedhadrons += p;
      }
    }
    if (chargedhadrons.empty()) return 0; //< leptonic tau
    if (pThadvis < 20*GeV) return 0; //< below threshold
    if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.55 : 1/50.;
    if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.40 : 1/110.;
    return 0;
  }


  /// ATLAS Run 1 tau smearing
  /// @todo Currently a copy of the crappy jet smearing that is probably wrong...
  inline Particle TAU_SMEAR_ATLAS_RUN1(const Particle& t) {
    // Const fractional resolution for now
    static const double resolution = 0.03;

    // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
    /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());
    normal_distribution<> d(1., resolution);
    const double fsmear = max(d(gen), 0.);
    const double mass = t.mass2() > 0 ? t.mass() : 0; //< numerical carefulness...
    return Particle(t.pid(), FourMomentum::mkXYZM(t.px()*fsmear, t.py()*fsmear, t.pz()*fsmear, mass));
  }


  /// ATLAS Run 2 tau smearing
  /// @todo Currently a copy of the Run 1 version
  inline Particle TAU_SMEAR_ATLAS_RUN2(const Particle& t) {
    return TAU_SMEAR_ATLAS_RUN1(t);
  }


  /// CMS Run 2 tau efficiency
  ///
  /// @todo Needs work; this is the dumb version from Delphes 3.3.2
  inline double TAU_EFF_CMS_RUN2(const Particle& t) {
    return (t.abspid() == PID::TAU) ? 0.6 : 0;
  }

  /// CMS Run 1 tau efficiency
  ///
  /// @todo Needs work; this is just a copy of the Run 2 version in Delphes 3.3.2
  inline double TAU_EFF_CMS_RUN1(const Particle& t) {
    return TAU_EFF_CMS_RUN2(t);
  }


  /// CMS Run 1 tau smearing
  /// @todo Currently a copy of the crappy ATLAS one
  inline Particle TAU_SMEAR_CMS_RUN1(const Particle& t) {
    return TAU_SMEAR_ATLAS_RUN1(t);
  }


  /// CMS Run 2 tau smearing
  /// @todo Currently a copy of the Run 1 version
  inline Particle TAU_SMEAR_CMS_RUN2(const Particle& t) {
    return TAU_SMEAR_CMS_RUN1(t);
  }

  //@}


  /// @name Jet efficiency and smearing functions
  //@{

  /// Return a constant 0 given a Jet as argument
  inline double JET_EFF_ZERO(const Jet& p) { return 0; }
  /// Return a constant 1 given a Jet as argument
  inline double JET_EFF_ONE(const Jet& p) { return 1; }

  /// Return 1 if the given Jet contains a b, otherwise 0
  inline double JET_BTAG_PERFECT(const Jet& j) { return j.bTagged() ? 1 : 0; }
  /// Return the ATLAS Run 1 jet flavour tagging efficiency for the given Jet
  inline double JET_BTAG_ATLAS_RUN1(const Jet& j) {
    if (j.bTagged()) return 0.80*tanh(0.003*j.pT()/GeV)*(30/(1+0.086*j.pT()/GeV));
    if (j.cTagged()) return 0.20*tanh(0.02*j.pT()/GeV)*(1/(1+0.0034*j.pT()/GeV));
    return 0.002 + 7.3e-6*j.pT()/GeV;
  }

  /// Return 1 if the given Jet contains a c, otherwise 0
  inline double JET_CTAG_PERFECT(const Jet& j) { return j.cTagged() ? 1 : 0; }

  /// Take a jet and return an unmodified copy
  /// @todo Modify constituent particle vectors for consistency
  /// @todo Set a null PseudoJet if the Jet is smeared?
  inline Jet JET_SMEAR_IDENTITY(const Jet& j) { return j; }

  /// ATLAS Run 1 jet smearing
  /// @todo This is a cluster-level flat 3% resolution, I think, and smearing is suboptimal: improve!
  inline Jet JET_SMEAR_ATLAS_RUN1(const Jet& j) {
    // Const fractional resolution for now
    static const double resolution = 0.03;

    // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
    /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
    /// @todo Need to isolate random generators to a single thread
    static random_device rd;
    static mt19937 gen(rd());
    normal_distribution<> d(1., resolution);
    const double fsmear = max(d(gen), 0.);
    const double mass = j.mass2() > 0 ? j.mass() : 0; //< numerical carefulness...
    return Jet(FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, mass));
  }

  /// CMS Run 1 jet smearing
  /// @todo Just a copy of the suboptimal ATLAS one: improve!!
  inline Jet JET_SMEAR_CMS_RUN1(const Jet& j) {
    return JET_SMEAR_ATLAS_RUN1(j);
  }

  //@}


}

#endif

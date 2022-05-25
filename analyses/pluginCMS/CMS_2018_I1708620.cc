/// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Forward energy flow at 13 TeV with CMS
  ///
  /// @note Rivet 3 conversion by A Buckley
  class CMS_2018_I1708620 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1708620);


    /// Initialise
    void init() {

      declare(FinalState(), "FS");
      const ChargedFinalState cfsBSCplus(Cuts::eta > 3.9 && Cuts::eta < 4.4);
      declare(cfsBSCplus, "cfsBSCplus");
      const ChargedFinalState cfsBSCminus(Cuts::eta > -4.4 && Cuts::eta < -3.9);
      declare(cfsBSCminus, "cfsBSCminus");

      book(_noe_inel, "TMP/Ninel");
      book(_noe_nsd, "TMP/Nnsd");
      book(_noe_bsc, "TMP/Nbsc");
      book(_noe_sd, "TMP/Nsd");
      book(_noe_nsd_sd, "TMP/Nnsdsd");

      book(_h_inel, 1, 1, 1);
      book(_h_nsd , 2, 1, 1);
      book(_h_et  , 3, 1, 1);
      book(_h_sd  , 4, 1, 1);
    }


    void analyze(const Event& event) {

      const ChargedFinalState& cfsBSCplus = applyProjection<ChargedFinalState>(event, "cfsBSCplus");
      const ChargedFinalState& cfsBSCminus = applyProjection<ChargedFinalState>(event, "cfsBSCminus");
      const bool bscplus = !cfsBSCplus.empty();
      const bool bscminus = !cfsBSCminus.empty();

      // Find final-state particles
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      // const Particles particlesByRapidity = fs.particlesByPt();
      // sortBy(particlesByRapidity, cmpMomByRap);
      const Particles particlesByRapidity = fs.particles(cmpMomByRap);
      const size_t num_particles = particlesByRapidity.size();

      // Find gaps, and choose the middle one as the event gap centre
      vector<double> gaps, midpoints;
      for (size_t ip = 1; ip < num_particles; ++ip) {
        const Particle& p1 = particlesByRapidity[ip-1];
        const Particle& p2 = particlesByRapidity[ip];
        const double gap = p2.rapidity() - p1.rapidity();
        const double mid = (p2.rapidity() + p1.rapidity()) / 2;
        gaps.push_back(gap);
        midpoints.push_back(mid);
      }
      const size_t imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
      const double gapcenter = midpoints[imid];

      // Assign particles to sides and compute Mx, My
      FourMomentum v4mx, v4my;
      for (const Particle& p : particlesByRapidity) {
        (p.rapidity() > gapcenter ? v4mx : v4my) += p.momentum();
      }
      const double mx = v4mx.mass();
      const double my = v4my.mass();

      // Compute xi variables
      const double xix = sqr(mx/(sqrtS()/GeV));
      const double xiy = sqr(my/(sqrtS()/GeV));
      const double xi_sd = max(xix, xiy);

      // Compute if inelastic, and other variables
      const bool inel = (xi_sd > 1e-6);
      if (inel) _noe_inel->fill();
      const bool bsc = (bscplus && bscminus);
      if (bsc) _noe_bsc->fill(); ///< @todo Not re-entry safe: FIX

      // Count/histogram backward and forward Et
      static const double YBEAM = 9.54;
      int nplus  = 0, nminus = 0;
      for (const Particle& p : particlesByRapidity) {
        const double eta = p.eta();
        if (!p.isVisible()) continue;
        if (inRange(eta, 2.866, 5.205)) nplus += 1;
        if (inRange(eta, -5.205, -2.866)) nminus += 1;
        if (bsc) _h_et->fill(eta-YBEAM, p.Et()/GeV);
      }

      // Categorise as NSD
      const bool nsd = (nminus > 0 && nplus > 0);
      if (nsd) _noe_nsd->fill();

      // Histogram
      for (const Particle& p : particlesByRapidity) {
        if (inel) _h_inel->fill(p.abseta(), p.E()/GeV);
        if (nsd) _h_nsd->fill(p.abseta(), p.E()/GeV);
      }

      // SD selection
      static const double ETAMIN = 3.152;
      static const double ETAMAX = 5.205;
      static const double EMIN   = 5*GeV;
      bool stableParticleEnergyCutMinus = false, stableParticleEnergyCutPlus = false;
      for (const Particle& p : particlesByRapidity) {
        if (p.E() < EMIN) continue;
        if (inRange(p.eta(), -ETAMAX, -ETAMIN)) stableParticleEnergyCutMinus = true;
        if (inRange(p.eta(),  ETAMIN,  ETAMAX)) stableParticleEnergyCutPlus = true;
      }

      // Select SD-enhanced events with the following condition
      const bool sd = (stableParticleEnergyCutPlus != stableParticleEnergyCutMinus); //< bool XOR
      if (sd) {
        _noe_sd->fill();
        for (const Particle& p : particlesByRapidity) {
          if (inRange(p.abseta(), ETAMIN, ETAMAX)) {
            if (stableParticleEnergyCutPlus && p.eta() > 0) _h_sd->fill(p.abseta(), p.E()/GeV);
            if (stableParticleEnergyCutMinus && p.eta() < 0) _h_sd->fill(p.abseta(), p.E()/GeV);
          } else { // CASTOR
            _h_sd->fill(p.abseta(), p.E()/GeV/2);
          }
        }
      }

      // Count how many are NSD and SD
      if (nsd && sd) _noe_nsd_sd->fill();
    }


    void finalize() {
      scale(_h_inel, 0.5/_noe_inel->sumW());
      scale(_h_nsd, 0.5/_noe_nsd->sumW());
      scale(_h_et, 1/_noe_bsc->sumW());
      scale(_h_sd, 1/_noe_sd->sumW());
      MSG_INFO( "Number of events of INEL: " << _noe_inel->effNumEntries() );
      MSG_INFO( "Number of events of NSD: " << _noe_nsd->effNumEntries() );
      MSG_INFO( "Number of events of SD: " << _noe_sd->effNumEntries() );
      MSG_INFO( "Number of events of NSD and SD contribution: " << _noe_nsd_sd->effNumEntries() );
    }


    Histo1DPtr _h_inel, _h_nsd, _h_et, _h_sd;
    CounterPtr _noe_inel, _noe_nsd, _noe_bsc, _noe_sd, _noe_nsd_sd;

  };


  RIVET_DECLARE_PLUGIN(CMS_2018_I1708620);

}

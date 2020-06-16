// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {

  


  /// @brief CDF Run II underlying event in Drell-Yan
  /// @author Hendrik Hoeth
  ///
  /// Measurement of the underlying event in Drell-Yan
  /// \f$ Z/\gamma^* \to e^+ e^- \f$ and
  /// \f$ Z/\gamma^* \to \mu^+ \mu^- \f$ events. The reconstructed
  /// Z defines the \f$ \phi \f$ orientation. A Z mass window cut is applied.
  ///
  /// @par Run conditions
  ///
  /// @arg \f$ \sqrt{s} = \f$ 1960 GeV
  /// @arg produce Drell-Yan events
  /// @arg Set particles with c*tau > 10 mm stable
  /// @arg Z decay mode: Z -> e+e- and Z -> mu+mu-
  /// @arg gamma decay mode: gamma -> e+e- and gamma -> mu+mu-
  /// @arg minimum invariant mass of the fermion pair coming from the Z/gamma: 70 GeV
  class CDF_2010_S8591881_DY : public Analysis {
  public:

    /// Constructor
    CDF_2010_S8591881_DY() : Analysis("CDF_2010_S8591881_DY")
    {
    }


    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      const ChargedFinalState cfs((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  0.5*GeV));
      const ChargedFinalState clfs((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  20*GeV));
      declare(cfs, "FS");
      declare(ChargedLeptons(clfs), "CL");

      // Book histograms
      book(_hist_tnchg      , 1, 1, 1);
      book(_hist_pnchg      , 1, 1, 2);
      book(_hist_anchg      , 1, 1, 3);
      book(_hist_pmaxnchg   , 2, 1, 1);
      book(_hist_pminnchg   , 2, 1, 2);
      book(_hist_pdifnchg   , 2, 1, 3);
      book(_hist_tcptsum    , 3, 1, 1);
      book(_hist_pcptsum    , 3, 1, 2);
      book(_hist_acptsum    , 3, 1, 3);
      book(_hist_pmaxcptsum , 4, 1, 1);
      book(_hist_pmincptsum , 4, 1, 2);
      book(_hist_pdifcptsum , 4, 1, 3);
      book(_hist_tcptave    , 5, 1, 1);
      book(_hist_pcptave    , 5, 1, 2);
      book(_hist_tcptmax    , 6, 1, 1);
      book(_hist_pcptmax    , 6, 1, 2);
      book(_hist_zptvsnchg  , 7, 1, 1);
      book(_hist_cptavevsnchg , 8, 1, 1);
      book(_hist_cptavevsnchgsmallzpt , 9, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {

      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 1) {
        MSG_DEBUG("Failed multiplicity cut");
        vetoEvent;
      }

      // Get the leptons
      const Particles& leptons = apply<ChargedLeptons>(e, "CL").chargedLeptons();

      // We want exactly two leptons of the same flavour.
      MSG_DEBUG("lepton multiplicity = " << leptons.size());
      if (leptons.size() != 2 || leptons[0].pid() != -leptons[1].pid() ) vetoEvent;

      // Lepton pT > 20 GeV
      if (leptons[0].pT()/GeV <= 20 || leptons[1].pT()/GeV <= 20) vetoEvent;

      // Lepton pair should have an invariant mass between 70 and 110 and |eta| < 6
      const FourMomentum dilepton = leptons[0].momentum() + leptons[1].momentum();
      if (!inRange(dilepton.mass()/GeV, 70., 110.) || fabs(dilepton.eta()) >= 6) vetoEvent;
      MSG_DEBUG("Dilepton mass = " << dilepton.mass()/GeV << " GeV");
      MSG_DEBUG("Dilepton pT   = " << dilepton.pT()/GeV << " GeV");

      // Calculate the observables
      size_t   numToward(0),     numAway(0);
      long int numTrans1(0),     numTrans2(0);
      double ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
      double ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);
      const double phiZ = dilepton.azimuthalAngle();
      const double pTZ  = dilepton.pT();
      /// @todo Replace with for
      for (Particles::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
        // Don't use the leptons
        /// @todo Replace with PID::isLepton
        if (abs(p->pid()) < 20) continue;

        const double dPhi = deltaPhi(p->momentum().phi(), phiZ);
        const double pT = p->pT();
        double rotatedphi = p->momentum().phi() - phiZ;
        while (rotatedphi < 0) rotatedphi += 2*PI;

        if (dPhi < PI/3.0) {
          ptSumToward += pT;
          ++numToward;
          if (pT > ptMaxToward)
            ptMaxToward = pT;
        } else if (dPhi < 2*PI/3.0) {
          if (rotatedphi <= PI) {
            ptSumTrans1 += pT;
            ++numTrans1;
            if (pT > ptMaxTrans1)
              ptMaxTrans1 = pT;
          }
          else {
            ptSumTrans2 += pT;
            ++numTrans2;
            if (pT > ptMaxTrans2)
              ptMaxTrans2 = pT;
          }
        } else {
          ptSumAway += pT;
          ++numAway;
          if (pT > ptMaxAway)
            ptMaxAway = pT;
        }
        // We need to subtract the two leptons from the number of particles to get the correct multiplicity
        _hist_cptavevsnchg->fill(numParticles-2, pT);
        if (pTZ < 10)
          _hist_cptavevsnchgsmallzpt->fill(numParticles-2, pT);
      }

      // Fill the histograms
      _hist_tnchg->fill(pTZ, numToward/(4*PI/3));
      _hist_pnchg->fill(pTZ, (numTrans1+numTrans2)/(4*PI/3));
      _hist_pmaxnchg->fill(pTZ, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _hist_pminnchg->fill(pTZ, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _hist_pdifnchg->fill(pTZ, abs(numTrans1-numTrans2)/(2*PI/3));
      _hist_anchg->fill(pTZ, numAway/(4*PI/3));

      _hist_tcptsum->fill(pTZ, ptSumToward/(4*PI/3));
      _hist_pcptsum->fill(pTZ, (ptSumTrans1+ptSumTrans2)/(4*PI/3));
      _hist_pmaxcptsum->fill(pTZ, (ptSumTrans1>ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3));
      _hist_pmincptsum->fill(pTZ, (ptSumTrans1<ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3));
      _hist_pdifcptsum->fill(pTZ, fabs(ptSumTrans1-ptSumTrans2)/(2*PI/3));
      _hist_acptsum->fill(pTZ, ptSumAway/(4*PI/3));

      if (numToward > 0) {
        _hist_tcptave->fill(pTZ, ptSumToward/numToward);
        _hist_tcptmax->fill(pTZ, ptMaxToward);
      }
      if ((numTrans1+numTrans2) > 0) {
        _hist_pcptave->fill(pTZ, (ptSumTrans1+ptSumTrans2)/(numTrans1+numTrans2));
        _hist_pcptmax->fill(pTZ, (ptMaxTrans1 > ptMaxTrans2 ? ptMaxTrans1 : ptMaxTrans2));
      }

      // We need to subtract the two leptons from the number of particles to get the correct multiplicity
      _hist_zptvsnchg->fill(numParticles-2, pTZ);
    }


    void finalize() {
    }

    //@}


  private:

    Profile1DPtr _hist_tnchg;
    Profile1DPtr _hist_pnchg;
    Profile1DPtr _hist_pmaxnchg;
    Profile1DPtr _hist_pminnchg;
    Profile1DPtr _hist_pdifnchg;
    Profile1DPtr _hist_anchg;
    Profile1DPtr _hist_tcptsum;
    Profile1DPtr _hist_pcptsum;
    Profile1DPtr _hist_pmaxcptsum;
    Profile1DPtr _hist_pmincptsum;
    Profile1DPtr _hist_pdifcptsum;
    Profile1DPtr _hist_acptsum;
    Profile1DPtr _hist_tcptave;
    Profile1DPtr _hist_pcptave;
    Profile1DPtr _hist_tcptmax;
    Profile1DPtr _hist_pcptmax;
    Profile1DPtr _hist_zptvsnchg;
    Profile1DPtr _hist_cptavevsnchg;
    Profile1DPtr _hist_cptavevsnchgsmallzpt;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2010_S8591881_DY);

}

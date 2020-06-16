// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// ATLAS jet substructure measurement
  class ATLAS_2012_I1094564 : public Analysis {
  public:

    ATLAS_2012_I1094564()
      : Analysis("ATLAS_2012_I1094564")
    {}


    // Returns constituents to make it easier to do the filtering
    PseudoJets splitjet(fastjet::PseudoJet jet, double& last_R, const FastJets& fj, bool& unclustered) const {

      // Build a new cluster sequence just using the constituents of this jet.
      fastjet::ClusterSequence cs(fj.clusterSeq()->constituents(jet), fastjet::JetDefinition(fastjet::cambridge_algorithm, M_PI/2.));

      // Get the jet back again
      vector<fastjet::PseudoJet> remadeJets = cs.inclusive_jets(0.);

      if ( remadeJets.size() != 1 ) return remadeJets;

      fastjet::PseudoJet remadeJet = remadeJets[0];
      fastjet::PseudoJet parent1, parent2;
      unclustered = false;

      while ( cs.has_parents(remadeJet, parent1, parent2) ) {
        if (parent1.squared_distance(parent2) < 0.09) break;
        if (parent1.m2() < parent2.m2()) {
          fastjet::PseudoJet tmp;
          tmp = parent1; parent1 = parent2; parent2 = tmp;
        }

        double ktdist = parent1.kt_distance(parent2);
        double rtycut2 = 0.3*0.3;
        if (parent1.m() < 0.67*remadeJet.m() && ktdist > rtycut2*remadeJet.m2()) {
          unclustered = true;
          break;
        } else {
          remadeJet = parent1;
        }
      }

      last_R = 0.5 * sqrt(parent1.squared_distance(parent2));
      return cs.constituents(remadeJet);
    }


    fastjet::PseudoJet filterjet(PseudoJets jets, double& stingy_R, const double def_R) const {
      if (stingy_R == 0.0) stingy_R = def_R;
      stingy_R = def_R < stingy_R ? def_R : stingy_R;
      fastjet::JetDefinition stingy_jet_def(fastjet::cambridge_algorithm, stingy_R);
      fastjet::ClusterSequence scs(jets, stingy_jet_def);
      vector<fastjet::PseudoJet> stingy_jets = sorted_by_pt(scs.inclusive_jets(0));
      fastjet::PseudoJet reconst_jet(0, 0, 0, 0);
      for (size_t isj = 0; isj < std::min((size_t) 3, stingy_jets.size()); ++isj) {
        reconst_jet += stingy_jets[isj];
      }
      return reconst_jet;
    }


    // These are custom functions for n-subjettiness.
    PseudoJets jetGetAxes(int n_jets, const PseudoJets& inputJets, double subR) const {
      // Sanity check
      if (inputJets.size() < (size_t) n_jets) {
        MSG_ERROR("Not enough input particles.");
        return inputJets;
      }

      // Get subjets, return
      fastjet::ClusterSequence sub_clust_seq(inputJets, fastjet::JetDefinition(fastjet::kt_algorithm, subR, fastjet::E_scheme, fastjet::Best));
      return sub_clust_seq.exclusive_jets(n_jets);
    }


    double jetTauValue(double beta, double jet_rad, const PseudoJets& particles, const PseudoJets& axes, double Rcut) const {
      double tauNum = 0.0;
      double tauDen = 0.0;

      if (particles.size() == 0) return 0.0;

      for (size_t i = 0; i < particles.size(); i++) {
        // find minimum distance (set R large to begin)
        double minR = 10000.0;
        for (size_t j = 0; j < axes.size(); j++) {
          double tempR = sqrt(particles[i].squared_distance(axes[j]));
          if (tempR < minR) minR = tempR;
        }
        if (minR > Rcut) minR = Rcut;
        // calculate nominator and denominator
        tauNum += particles[i].perp() * pow(minR,beta);
        tauDen += particles[i].perp() * pow(jet_rad,beta);
      }

      // return N-subjettiness (or 0 if denominator is 0)
      return safediv(tauNum, tauDen, 0);
    }


    void jetUpdateAxes(double beta, const PseudoJets& particles, PseudoJets& axes) const {
      vector<int> belongsto;
      for (size_t i = 0; i < particles.size(); i++) {
        // find minimum distance axis
        int assign = 0;
        double minR = 10000.0;
        for (size_t j = 0; j < axes.size(); j++) {
          double tempR = sqrt(particles[i].squared_distance(axes[j]));
          if (tempR < minR) {
            minR = tempR;
            assign = j;
          }
        }
        belongsto.push_back(assign);
      }

      // iterative step
      double deltaR2, distphi;
      vector<double> ynom, phinom, den;

      ynom.resize(axes.size());
      phinom.resize(axes.size());
      den.resize(axes.size());

      for (size_t i = 0; i < particles.size(); i++) {
        distphi = particles[i].phi() - axes[belongsto[i]].phi();
        deltaR2 = particles[i].squared_distance(axes[belongsto[i]]);
        if (deltaR2 == 0.) continue;
        if (abs(distphi) <= M_PI) phinom.at(belongsto[i]) += particles[i].perp() * particles[i].phi() * pow(deltaR2, (beta-2)/2);
        else if ( distphi > M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (-2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);
        else if ( distphi < M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (+2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);
        ynom.at(belongsto[i]) += particles[i].perp() * particles[i].rap() * pow(deltaR2, (beta-2)/2);
        den.at(belongsto[i])  += particles[i].perp() * pow(deltaR2, (beta-2)/2);
      }

      // reset to new axes
      for (size_t j = 0; j < axes.size(); j++) {
        if (den[j] == 0.) axes.at(j) = axes[j];
        else {
          double phi_new = fmod( 2*M_PI + (phinom[j] / den[j]), 2*M_PI );
          double pt_new  = axes[j].perp();
          double y_new   = ynom[j] / den[j];
          double px = pt_new * cos(phi_new);
          double py = pt_new * sin(phi_new);
          double pz = pt_new * sinh(y_new);
          axes.at(j).reset(px, py, pz, axes[j].perp()/2);
        }
      }
    }


    void init() {

      /// Projections:
      FinalState fs((Cuts::etaIn(-4.5, 4.5) && Cuts::pT >=  0.*GeV));
      declare(fs, "FS");
      declare(FastJets(fs, FastJets::ANTIKT, 1.0), "AKT");
      declare(FastJets(fs, FastJets::CAM, 1.2)   , "CA" );

      /// Histograms:
      {Histo1DPtr tmp; _h_camass.add(200, 300, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_camass.add(300, 400, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_camass.add(400, 500, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_camass.add(500, 600, book(tmp, 4, 1, 1));}

      {Histo1DPtr tmp; _h_filtmass.add(200, 300, book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_filtmass.add(300, 400, book(tmp, 6, 1, 1));}
      {Histo1DPtr tmp; _h_filtmass.add(400, 500, book(tmp, 7, 1, 1));}
      {Histo1DPtr tmp; _h_filtmass.add(500, 600, book(tmp, 8, 1, 1));}

      {Histo1DPtr tmp; _h_ktmass.add(200, 300, book(tmp,  9, 1, 1));}
      {Histo1DPtr tmp; _h_ktmass.add(300, 400, book(tmp, 10, 1, 1));}
      {Histo1DPtr tmp; _h_ktmass.add(400, 500, book(tmp, 11, 1, 1));}
      {Histo1DPtr tmp; _h_ktmass.add(500, 600, book(tmp, 12, 1, 1));}

      {Histo1DPtr tmp; _h_ktd12.add(200, 300, book(tmp, 13, 1, 1));}
      {Histo1DPtr tmp; _h_ktd12.add(300, 400, book(tmp, 14, 1, 1));}
      {Histo1DPtr tmp; _h_ktd12.add(400, 500, book(tmp, 15, 1, 1));}
      {Histo1DPtr tmp; _h_ktd12.add(500, 600, book(tmp, 16, 1, 1));}

      {Histo1DPtr tmp; _h_ktd23.add(200, 300, book(tmp, 17, 1 ,1));}
      {Histo1DPtr tmp; _h_ktd23.add(300, 400, book(tmp, 18, 1 ,1));}
      {Histo1DPtr tmp; _h_ktd23.add(400, 500, book(tmp, 19, 1 ,1));}
      {Histo1DPtr tmp; _h_ktd23.add(500, 600, book(tmp, 20, 1 ,1));}

      {Histo1DPtr tmp; _h_cat21.add(200, 300, book(tmp, 21, 1, 1));}
      {Histo1DPtr tmp; _h_cat21.add(300, 400, book(tmp, 22, 1, 1));}
      {Histo1DPtr tmp; _h_cat21.add(400, 500, book(tmp, 23, 1, 1));}
      {Histo1DPtr tmp; _h_cat21.add(500, 600, book(tmp, 24, 1, 1));}

      {Histo1DPtr tmp; _h_cat32.add(200, 300, book(tmp, 25, 1, 1));}
      {Histo1DPtr tmp; _h_cat32.add(300, 400, book(tmp, 26, 1, 1));}
      {Histo1DPtr tmp; _h_cat32.add(400, 500, book(tmp, 27, 1, 1));}
      {Histo1DPtr tmp; _h_cat32.add(500, 600, book(tmp, 28, 1, 1));}

      {Histo1DPtr tmp; _h_ktt21.add(200, 300, book(tmp, 29, 1, 1));}
      {Histo1DPtr tmp; _h_ktt21.add(300, 400, book(tmp, 30, 1, 1));}
      {Histo1DPtr tmp; _h_ktt21.add(400, 500, book(tmp, 31, 1, 1));}
      {Histo1DPtr tmp; _h_ktt21.add(500, 600, book(tmp, 32, 1, 1));}

      {Histo1DPtr tmp; _h_ktt32.add(200, 300, book(tmp, 33, 1, 1));}
      {Histo1DPtr tmp; _h_ktt32.add(300, 400, book(tmp, 34, 1, 1));}
      {Histo1DPtr tmp; _h_ktt32.add(400, 500, book(tmp, 35, 1, 1));}
      {Histo1DPtr tmp; _h_ktt32.add(500, 600, book(tmp, 36, 1, 1));}
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      using namespace fastjet;

      // Get anti-kt jets with p_T > 200 GeV, check abs(y) < 2, and fill mass histograms
      const FastJets& ktfj = apply<FastJets>(event, "AKT");
      PseudoJets ktjets = ktfj.pseudoJetsByPt(200*GeV);
      for (const PseudoJet ajet : ktjets) {
        if (abs(ajet.rap()) < 2) {
          _h_ktmass.fill(ajet.perp(), ajet.m(), weight);
        }
      }

      // Same as above but C/A jets
      const FastJets& cafj = apply<FastJets>(event, "CA");
      PseudoJets cajets = cafj.pseudoJetsByPt(200*GeV);
      for (const PseudoJet ajet : cajets) {
        if (abs(ajet.rap()) < 2) {
          _h_camass.fill(ajet.perp(), ajet.m(), weight);
        }
      }

      // Split and filter.
      // Only do this to C/A jets in this analysis.
      for (const PseudoJet pjet : cajets) {
        if ( pjet.perp() > 600 || abs(pjet.rap()) > 2) continue;
        double dR = 0;
        bool unclustered = false;
        PseudoJets split_jets = splitjet(pjet, dR, cafj, unclustered);
        if ( (dR < 0.15) || (unclustered == false) ) continue;
        PseudoJet filt_jet = filterjet(split_jets, dR, 0.3);
        _h_filtmass.fill(filt_jet.perp(), filt_jet.m(), weight);
      }

      // Use the two last stages of clustering to get sqrt(d_12) and sqrt(d_23).
      // Only use anti-kt jets in this analysis.
      for (const PseudoJet pjet : ktjets) {
        if (pjet.perp() > 600 || abs(pjet.rap()) > 2) continue;
        ClusterSequence subjet_cseq(ktfj.clusterSeq()->constituents(pjet), JetDefinition(kt_algorithm, M_PI/2.));
        double d_12 = subjet_cseq.exclusive_dmerge(1) * M_PI*M_PI/4.;
        double d_23 = subjet_cseq.exclusive_dmerge(2) * M_PI*M_PI/4.;
        _h_ktd12.fill(pjet.perp(), sqrt(d_12), weight);
        _h_ktd23.fill(pjet.perp(), sqrt(d_23), weight);
      }

      // N-subjettiness, use beta = 1 (no rationale given).
      // Uses the functions defined above.
      // C/A jets first, anti-kt after.
      double beta = 1.;
      //Rcut is used for particles that are very far from the closest axis. At 10
      //is has no impact on the outcome of the calculation
      double Rcut = 10.;
      for (const PseudoJet pjet : cajets) {
        if (pjet.perp() > 600*GeV || fabs(pjet.rap()) > 2) continue;

        const PseudoJets constituents = cafj.clusterSeq()->constituents(pjet);
        if (constituents.size() < 3) continue;

        const PseudoJets axis1 = jetGetAxes(1, constituents, M_PI/2.0);
        const PseudoJets axis2 = jetGetAxes(2, constituents, M_PI/2.0);
        const PseudoJets axis3 = jetGetAxes(3, constituents, M_PI/2.0);

        const double radius = 1.2;
        const double tau1 = jetTauValue(beta, radius, constituents, axis1, Rcut);
        const double tau2 = jetTauValue(beta, radius, constituents, axis2, Rcut);
        const double tau3 = jetTauValue(beta, radius, constituents, axis3, Rcut);

        if (tau1 == 0 || tau2 == 0) continue;
        _h_cat21.fill(pjet.perp(), tau2/tau1, weight);
        _h_cat32.fill(pjet.perp(), tau3/tau2, weight);
      }

      for (const PseudoJet pjet : ktjets) {
        if (pjet.perp() > 600*GeV || fabs(pjet.rap()) > 2) continue;

        const PseudoJets constituents = ktfj.clusterSeq()->constituents(pjet);
        if (constituents.size() < 3) continue;

        const PseudoJets axis1 = jetGetAxes(1, constituents, M_PI/2.0);
        const PseudoJets axis2 = jetGetAxes(2, constituents, M_PI/2.0);
        const PseudoJets axis3 = jetGetAxes(3, constituents, M_PI/2.0);

        const double radius = 1.0;
        const double tau1 = jetTauValue(beta, radius, constituents, axis1, Rcut);
        const double tau2 = jetTauValue(beta, radius, constituents, axis2, Rcut);
        const double tau3 = jetTauValue(beta, radius, constituents, axis3, Rcut);
        if (tau1 == 0 || tau2 == 0) continue;

        _h_ktt21.fill(pjet.perp(), tau2/tau1, weight);
        _h_ktt32.fill(pjet.perp(), tau3/tau2, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (Histo1DPtr h : _h_camass.histos()) normalize(h);
      for (Histo1DPtr h : _h_filtmass.histos()) normalize(h);
      for (Histo1DPtr h : _h_ktmass.histos()) normalize(h);
      for (Histo1DPtr h : _h_ktd12.histos()) normalize(h);
      for (Histo1DPtr h : _h_ktd23.histos()) normalize(h);
      for (Histo1DPtr h : _h_cat21.histos()) normalize(h);
      for (Histo1DPtr h : _h_cat32.histos()) normalize(h);
      for (Histo1DPtr h : _h_ktt21.histos()) normalize(h);
      for (Histo1DPtr h : _h_ktt32.histos()) normalize(h);
    }

  private:

    BinnedHistogram _h_camass;
    BinnedHistogram _h_filtmass;
    BinnedHistogram _h_ktmass;
    BinnedHistogram _h_ktd12;
    BinnedHistogram _h_ktd23;
    BinnedHistogram _h_cat21;
    BinnedHistogram _h_cat32;
    BinnedHistogram _h_ktt21;
    BinnedHistogram _h_ktt32;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1094564);

}

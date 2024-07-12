// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/EventMixingFinalState.hh"

namespace Rivet {


  /// @brief Angular correlations of identified particles in pp at 7 TeV.
  ///
  /// Also showcasing use of EventMixingFinalState.
  class ALICE_2016_I1507157 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2016_I1507157);


    /// @name Analysis methods
    /// @{

    /// @brief Calculate angular distance between particles.
    double phaseDif(double a1, double a2, const pair<double, double>& edges) {
      double dif = a1 - a2;
      while (dif < edges.first)
        dif += 2*M_PI;
      while (dif > edges.second)
        dif -= 2*M_PI;
      return dif;
    }


    /// @brief Get the minimal and maximal x values from a Scatter2D (refdata).
    pair<double, double> xEdges(const YODA::Scatter2D& h) {
      double xMin = 0;
      double xMax = 0;
      for (const auto& p : h.points()) {
        if (xMin > p.xMin()) xMin = p.xMin();
        if (xMax < p.xMax()) xMax = p.xMax();
      }
      return {xMin, xMax};
    }


    /// Book histograms and initialise projections before the run
    void init() {

      const double etamax = 0.8;
      const double pTmin = 0.2; // GeV
      const double pTmax = 2.5; //GeV

      // Trigger projection.
      declare(ALICE::V0AndTrigger(), "V0-AND");
      // Charged tracks used to manage the mixing observable.
      ChargedFinalState cfsMult(Cuts::abseta < etamax);
      declare(cfsMult, "CFSMult");

      // Primary particles.
      PrimaryParticles pp({Rivet::PID::PIPLUS, Rivet::PID::KPLUS,
	      Rivet::PID::K0S, Rivet::PID::K0L, Rivet::PID::PROTON,
	      Rivet::PID::NEUTRON, Rivet::PID::LAMBDA, Rivet::PID::SIGMAMINUS,
       	Rivet::PID::SIGMAPLUS, Rivet::PID::XIMINUS, Rivet::PID::XI0,
	      Rivet::PID::OMEGAMINUS},Cuts::abseta < etamax && 
        Cuts::pT > pTmin*GeV && Cuts::pT < pTmax*GeV);
      declare(pp,"APRIM");

      // The event mixing projection
      declare(EventMixingFinalState(cfsMult, pp, 5, 0, 100, 10, defaultWeightIndex()),"EVM");
      // The particle pairs.
      pid = {{211, -211}, {321, -321}, {2212, -2212}, {3122, -3122}, {211, 211},
             {321, 321}, {2212, 2212}, {3122, 3122}, {2212, 3122}, {2212, -3122}};
      // The differing pT cuts per pair, in GeV.
      pTcuts = {{0.2, 0.2},{0.3, 0.3},{0.5,0.5},{0.6,0.6},{0.2,0.2},
	      {0.3,0.3},{0.5,0.5},{0.6,0.6},{0.5,0.6},{0.5,0.6}};
      // The associated histograms in the data file.
      vector<string> refdata = {"d04-x01-y01","d04-x01-y02","d04-x01-y03",
        "d06-x01-y02","d05-x01-y01","d05-x01-y02","d05-x01-y03","d06-x01-y01",
        "d01-x01-y02","d02-x01-y02"};
      // Resize all the analysis object containers to right size.
      ratio.resize(refdata.size());
      signal.resize(refdata.size());
      background.resize(refdata.size());
      nsp.resize(refdata.size());
      nmp.resize(refdata.size());
      for (int i = 0, N = refdata.size(); i < N; ++i) {
        const YODA::Scatter2D& tmp = refData(refdata[i]);
        // The ratio plots.
        book(ratio[i], refdata[i], true);
        // Signal and mixed background should not be displayed.
        book(signal[i], "TMP/" + refdata[i] + "-s", tmp);
        book(background[i], "TMP/" + refdata[i] + "-b", tmp);
        // Number of signal and mixed pairs for normalization.
        book(nsp[i],"TMP/nsp"+std::to_string(i));
        book(nmp[i],"TMP/nmp"+std::to_string(i));
        // The differing deltaphi histogram edges per pair.
        deltaphi.push_back(xEdges(tmp));
      }
    }


    void fillPair(const Particle& p1, const Particle& p2, vector<Histo1DPtr>& histos, 
      vector<CounterPtr>& sow) {
	   if (isSame(p1,p2)) return;
          // If the pair is not within eta acceptance, we can continue early.
          if (abs(p1.eta() - p2.eta()) > 1.3) return;
          // Figure out which pid pair we are looking at.
          int iPair = -1;
          for (int i = 0, N = pid.size(); i < N; ++i) {
            if (pid[i].first == p1.pid() && pid[i].second == p2.pid()) {
              iPair = i;
              break;
            }
          }
          // If the pair is not in the analysis, don't fill anything.
          if (iPair < 0) return;
          // Apply min pT cuts, varies for different species.
          if (p1.pT() < pTcuts[iPair].first || p2.pT() < pTcuts[iPair].second) return;
          const double dPhi = phaseDif(p1.phi(), p2.phi(), deltaphi[iPair]);
          histos[iPair]->fill(dPhi);
          sow[iPair]->fill();
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Triggering.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")()) return;

      // The projections for signal and mixed event background.
      const PrimaryParticles& pp = 
        applyProjection<PrimaryParticles>(event,"APRIM");
      const EventMixingFinalState& evm = 
        applyProjection<EventMixingFinalState>(event, "EVM");

      // Test if we have enough mixing events available to continue.
      if (!evm.hasMixingEvents()) return;

      for (const Particle& p1 : pp.particles()) {
	      // First do the signal histograms. 
        for (const Particle& p2 : pp.particles())
	        fillPair(p1, p2, signal, nsp);
	      // Then do the background
        for (const Particle& p2 : evm.particles())
	        fillPair(p1, p2, background, nmp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int i = 0, N = pid.size(); i < N; ++i) {
	    // Scaling factor eqns. (2)-(5) in the paper.
        double sc = nmp[i]->sumW() / nsp[i]->sumW();
        signal[i]->scaleW(sc);
        divide(signal[i],background[i],ratio[i]);
      }
    }

    /// @}


    /// Analysis variables.
    vector<pair<int, int> > pid;
    vector<pair<double, double> > pTcuts;
    vector<pair<double, double> > deltaphi;
    /// @name Histograms and counters
    /// @{
    vector<Histo1DPtr> signal;
    vector<Histo1DPtr> background;
    vector<Scatter2DPtr> ratio;
    vector<CounterPtr> nsp;
    vector<CounterPtr> nmp;
    /// @}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2016_I1507157);

}

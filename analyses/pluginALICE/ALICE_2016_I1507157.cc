// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/EventMixingFinalState.hh"
namespace Rivet {


  /// @brief Correlations of identified particles in pp.
  /// Also showcasing use of EventMixingFinalState.
  class ALICE_2016_I1507157 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2016_I1507157);


    /// @name Analysis methods
    //@{

    /// @brief Calculate angular distance between particles.
    double phaseDif(double a1, double a2){
      double dif = a1 - a2;
      while (dif < -M_PI/2)
        dif += 2*M_PI;
      while (dif > 3*M_PI/2)
        dif -= 2*M_PI;
      return dif;
    }


    /// Book histograms and initialise projections before the run
    void init() {

      double etamax = 0.8;
      double pTmin = 0.5; // GeV

      // Trigger
      declare(ALICE::V0AndTrigger(), "V0-AND");
      // Charged tracks used to manage the mixing observable.
      ChargedFinalState cfsMult(Cuts::abseta < etamax);
      declare(cfsMult, "CFSMult");

      // Primary particles.
      PrimaryParticles pp({Rivet::PID::PIPLUS, Rivet::PID::KPLUS,
	Rivet::PID::K0S, Rivet::PID::K0L, Rivet::PID::PROTON,
	Rivet::PID::NEUTRON, Rivet::PID::LAMBDA, Rivet::PID::SIGMAMINUS,
       	Rivet::PID::SIGMAPLUS, Rivet::PID::XIMINUS, Rivet::PID::XI0,
	Rivet::PID::OMEGAMINUS},Cuts::abseta < etamax && Cuts::pT > pTmin*GeV);
      declare(pp,"APRIM");

      // The event mixing projection
      declare(EventMixingFinalState(cfsMult, pp, 5, 0, 100, 10, defaultWeightIndex()),"EVM");
      // The particle pairs.
      pid = {{211, -211}, {321, -321}, {2212, -2212}, {3122, -3122}, {211, 211},
             {321, 321}, {2212, 2212}, {3122, 3122}, {2212, 3122}, {2212, -3122}};
      // The associated histograms in the data file.
      vector<string> refdata = {"d04-x01-y01","d04-x01-y02","d04-x01-y03",
        "d06-x01-y02","d05-x01-y01","d05-x01-y02","d05-x01-y03","d06-x01-y01",
        "d01-x01-y02","d02-x01-y02"};
      ratio.resize(refdata.size());
      signal.resize(refdata.size());
      background.resize(refdata.size());
      for (int i = 0, N = refdata.size(); i < N; ++i) {
        // The ratio plots.
        book(ratio[i], refdata[i], true);
        // Signal and mixed background.
        book(signal[i], "/TMP/" + refdata[i] + "-s", refData(refdata[i]));
        book(background[i], "/TMP/" + refdata[i] + "-b", refData(refdata[i]));
        // Number of signal and mixed pairs.
        nsp.push_back(0.);
        nmp.push_back(0.);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Triggering
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")()) return;
      // The projections
      const PrimaryParticles& pp = 
        applyProjection<PrimaryParticles>(event,"APRIM");
      const EventMixingFinalState& evm = 
        applyProjection<EventMixingFinalState>(event, "EVM");
      // Test if we have enough mixing events available to continue.
      if (!evm.hasMixingEvents()) return;
      for(const Particle& p1 : pp.particles()) {
        // Start by doing the signal distributions
	for(const Particle& p2 : pp.particles()) {
	  if(isSame(p1,p2))
	    continue;
	  double dEta = abs(p1.eta() - p2.eta());
	  double dPhi = phaseDif(p1.phi(), p2.phi());
	  if(dEta < 1.3) {
	    for (int i = 0, N = pid.size(); i < N; ++i) {
	      int pid1 = pid[i].first;
	      int pid2 = pid[i].second;
	      bool samesign = (pid1 * pid2 > 0);
	      if (samesign && ((pid1 == p1.pid() && pid2 == p2.pid()) ||
		 (pid1 == -p1.pid() && pid2 == -p2.pid()))) {
	        signal[i]->fill(dPhi);
		nsp[i] += 1.0;
	      }
	      if (!samesign && abs(pid1) == abs(pid2) &&
		  pid1 == p1.pid() && pid2 == p2.pid()) {
	            signal[i]->fill(dPhi);
		    nsp[i] += 1.0;
	      }
	      if (!samesign && abs(pid1) != abs(pid2) &&
		  ( (pid1 == p1.pid() && pid2 == p2.pid()) ||
		  (pid2 == p1.pid() && pid1 == p2.pid()) ) ) {
	            signal[i]->fill(dPhi);
		    nsp[i] += 1.0;
	      }
	    }
	  }
	}
	// Then do the background distribution
	for(const Particle& pMix : evm.particles()){
	  double dEta = abs(p1.eta() - pMix.eta());
	  double dPhi = phaseDif(p1.phi(), pMix.phi());
	  if(dEta < 1.3) {
	    for (int i = 0, N = pid.size(); i < N; ++i) {
	      int pid1 = pid[i].first;
	      int pid2 = pid[i].second;
	      bool samesign = (pid1 * pid2 > 0);
	      if (samesign && ((pid1 == p1.pid() && pid2 == pMix.pid()) ||
		 (pid1 == -p1.pid() && pid2 == -pMix.pid()))) {
	            background[i]->fill(dPhi);
		    nmp[i] += 1.0;
	      }
	      if (!samesign && abs(pid1) == abs(pid2) &&
		  pid1 == p1.pid() && pid2 == pMix.pid()) {
	            background[i]->fill(dPhi);
		    nmp[i] += 1.0;
	      }
	      if (!samesign && abs(pid1) != abs(pid2) &&
		  ( (pid1 == p1.pid() && pid2 == pMix.pid()) ||
		  (pid2 == p1.pid() && pid1 == pMix.pid()) ) ) {
	            background[i]->fill(dPhi);
		    nmp[i] += 1.0;
	      }
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int i = 0, N = pid.size(); i < N; ++i) {
        double sc = nmp[i] / nsp[i];
        signal[i]->scaleW(sc);
        divide(signal[i],background[i],ratio[i]);
      }
    }

    //@}


    /// @name Histograms
    //@{
    vector<pair<int, int> > pid;
    vector<Histo1DPtr> signal;
    vector<Histo1DPtr> background;
    vector<Scatter2DPtr> ratio;
    vector<double> nsp;
    vector<double> nmp;

    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2016_I1507157);


}

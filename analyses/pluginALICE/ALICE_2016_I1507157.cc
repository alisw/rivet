// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/EventMixingFinalState.hh"
namespace Rivet {


  /// @brief ALICE correlations of identified particles in pp
  /// Also showcasing use of EventMixingFinalState.
  class ALICE_2016_I1507157 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2016_I1507157);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      double etamax = 0.8;
      double pTmin = 0.5; // GeV

      // Trigger
      declare(ALICE::V0AndTrigger(), "V0-AND");
      // Charged tracks used to manage the mixing observable.
      ChargedFinalState cfsMult(Cuts::abseta < etamax);
      addProjection(cfsMult, "CFSMult");

      // Primary particles.
      PrimaryParticles pp({Rivet::PID::PIPLUS, Rivet::PID::KPLUS,
            Rivet::PID::K0S, Rivet::PID::K0L, Rivet::PID::PROTON,
            Rivet::PID::NEUTRON, Rivet::PID::LAMBDA, Rivet::PID::SIGMAMINUS,
            Rivet::PID::SIGMAPLUS, Rivet::PID::XIMINUS, Rivet::PID::XI0,
            Rivet::PID::OMEGAMINUS},Cuts::abseta < etamax && Cuts::pT > pTmin*GeV);
      addProjection(pp,"APRIM");

      // The event mixing projection
      declare(EventMixingFinalState(cfsMult, pp, 5, 0, 100, 10),"EVM");
      // The particle pairs.
      pid = {{211, -211}, {321, -321}, {2212, -2212}, {3122, -3122}, {211, 211},
             {321, 321}, {2212, 2212}, {3122, 3122}, {2212, 3122}, {2212, -3122}};
      // The associated histograms in the data file.
      vector<string> refdata = {"d04-x01-y01","d04-x01-y02","d04-x01-y03",
                                "d06-x01-y02","d05-x01-y01","d05-x01-y02","d05-x01-y03","d06-x01-y01",
                                "d01-x01-y02","d02-x01-y02"};
      for (int i = 0, N = refdata.size(); i < N; ++i) {
        // The ratio plots.
        ratio.push_back(bookScatter2D(refdata[i], true));
        // Signal and mixed background.
        signal.push_back(bookHisto1D("/TMP/" + refdata[i] +
                                     "-s", *ratio[i], refdata[i] + "-s"));
        background.push_back(bookHisto1D("/TMP/" + refdata[i] +
                                         "-b", *ratio[i], refdata[i] + "-b"));
        // Number of signal and mixed pairs.
        nsp.push_back(0.);
        nmp.push_back(0.);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Triggering
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")()) return;

      // The primary particles
      const PrimaryParticles& pp = apply<PrimaryParticles>(event, "APRIM");
      const Particles pparticles = pp.particles();

      // The mixed events
      const EventMixingFinalState& evm = apply<EventMixingFinalState>(event, "EVM");
      const vector<Particles> mixEvents = evm.getMixingEvents();
      if (mixEvents.empty()) vetoEvent;

      // Make a vector of mixed event particles
      vector<Particle> mixParticles;
      size_t pSize = 0;
      for (size_t i = 0; i < mixEvents.size(); ++i)
        pSize += mixEvents[i].size();
      mixParticles.reserve(pSize);
      for (size_t i = 0; i < mixEvents.size(); ++i)
        mixParticles.insert(mixParticles.end(), mixEvents[i].begin(), mixEvents[i].end());
      random_shuffle(mixParticles.begin(), mixParticles.end());

      for (size_t ip1 = 0; ip1 < pparticles.size()-1; ++ip1) {
        const Particle& p1 = pparticles[ip1];

        // Start by doing the signal distributions
        for (size_t ip2 = 0; ip2 < pparticles.size(); ++ip1) {
          if (ip1 == ip2) continue;
          const Particle& p2 = pparticles[ip2];
          const double dEta = deltaEta(p1, p2);
          const double dPhi = deltaPhi(p1, p2, true);
          if (dEta > 1.3) continue;
          for (int i = 0, N = pid.size(); i < N; ++i) {
            const int pid1 = pid[i].first;
            const int pid2 = pid[i].second;
            const bool samesign = (pid1 * pid2 > 0);
            const bool pidmatch1 = (pid1 == p1.pid() && pid2 == p2.pid()) || (pid1 == -p1.pid() && pid2 == -p2.pid());
            const bool pidmatch2 = abs(pid1) == abs(pid2) && pid1 == p1.pid() && pid2 == p2.pid();
            const bool pidmatch3 = abs(pid1) != abs(pid2) && ( (pid1 == p1.pid() && pid2 == p2.pid()) || (pid2 == p1.pid() && pid1 == p2.pid()) );
            if ((samesign && pidmatch1) || (!samesign && (pidmatch2 || pidmatch3))) {
              signal[i]->fill(dPhi, weight);
              nsp[i] += 1.0;
            }
          }
        }

        // Then do the background distribution
        for (const Particle& pMix : mixParticles){
          const double dEta = deltaEta(p1, pMix);
          const double dPhi = deltaPhi(p1, pMix, true);
          if (dEta > 1.3) continue;
          for (int i = 0, N = pid.size(); i < N; ++i) {
            const int pid1 = pid[i].first;
            const int pid2 = pid[i].second;
            const bool samesign = (pid1 * pid2 > 0);
            const bool pidmatch1 = (pid1 == p1.pid() && pid2 == pMix.pid()) || (pid1 == -p1.pid() && pid2 == -pMix.pid());
            const bool pidmatch2 = abs(pid1) == abs(pid2) && pid1 == p1.pid() && pid2 == pMix.pid();
            const bool pidmatch3 = abs(pid1) != abs(pid2) && ( (pid1 == p1.pid() && pid2 == pMix.pid()) || (pid2 == p1.pid() && pid1 == pMix.pid()) );
            if ((samesign && pidmatch1) || (!samesign && (pidmatch2 || pidmatch3))) {
              background[i]->fill(dPhi, weight);
              nmp[i] += 1.0;
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int i = 0, N = pid.size(); i < N; ++i) {
        const double sc = nmp[i] / nsp[i];
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


  DECLARE_RIVET_PLUGIN(ALICE_2016_I1507157);

}

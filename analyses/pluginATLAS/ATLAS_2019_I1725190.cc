// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Smearing.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {


  /// @brief Dilepton high-mass resonance search
  ///
  /// @todo Use the proper smearing system rather than hand-rolled sampling
  class ATLAS_2019_I1725190 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2019_I1725190);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON && Cuts::Et > 30*GeV &&
                                 (Cuts::abseta < 2.47 && !Cuts::absetaIn(1.37, 1.52)));
      SmearedParticles recoelectrons(electrons, PARTICLE_EFF_ONE); //ELECTRON_EFF_ATLAS_RUN2_MEDIUM);
      declare(recoelectrons, "Elecs");

      PromptFinalState muons(Cuts::abspid == PID::MUON && Cuts::pT > 30*GeV && Cuts::abseta < 2.5);
      SmearedParticles recomuons(muons, PARTICLE_EFF_ONE);
      // [](const Particle& m) -> double {
      //   if (m.pT() < 1*TeV) return 0.69;
      //   if (m.pT() > 2.5*TeV) return 0.57;
      //   return 0.69 - 0.12*(m.pT() - 1*TeV)/(2.5*TeV - 1*TeV);
      // });
      declare(recomuons, "Muons");

      // Book histograms
      book(_h_mee, 1, 1, 1);
      book(_h_mmm, 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get leptons
      const Particles elecs = apply<ParticleFinder>(event, "Elecs").particlesByPt();
      const Particles muons = apply<ParticleFinder>(event, "Muons").particlesByPt();

      // MSG_INFO("Num e, mu = " << elecs.size() << ", " << muons.size());
      // for (const Particle& e : elecs) MSG_INFO(e.mom());
      // for (const Particle& m : muons) MSG_INFO(m.mom());

      // Isolation
      /// @todo Can't be done from provided information?
      // Particles isoelecs, isomuons;
      // for (const Particle& e : elecs) {
      // }

      // Choose the highest-pT lepton pair, preferring electrons
      if (elecs.size() < 2 && muons.size() < 2) vetoEvent;
      const Particle l1 = (elecs.size() >= 2) ? elecs[0] : muons[0];
      const Particle l2 = (elecs.size() >= 2) ? elecs[1] : muons[1];

      // Require opposite sign for muons only
      const bool mumu = (l1.abspid() == PID::MUON);
      if (mumu && l1.pid()*l2.pid() > 0) vetoEvent;

      // Get the true dilepton pair
      const FourMomentum pll = l1.mom() + l2.mom();
      const double mll = pll.mass();

      // Make sure we're in a region where the smearing and efficiencies are well-behaved
      if (mll < 200*GeV) vetoEvent;

      // Apply dilepton efficiency curves (as function of mX ~ mll)
      const double eff = mumu ?
        (0.54 - (mll - 200*GeV)/(6000*GeV - 200*GeV) * (0.54 - 0.38))  :
        (0.74 - 0.04*exp(-(mll-200*GeV)/100*GeV) - 0.08*exp(-(mll-200*GeV)/1000*GeV));
      if (rand01() > eff) vetoEvent;

      // Smear the dilepton mass with a CB + Gauss function
      double muCB, sigCB, alpCB, nCB, muG, sigG, kappa;
      if (!mumu) {
        const double lnm = log(mll);
        static const vector<double> pmuCB = {0.13287, -0.410663, -0.0126743, 2.9547e-6};
        muCB = pmuCB[0] + pmuCB[1]/lnm + pmuCB[2]*lnm + pmuCB[3]*pow(lnm, 4);
        static const vector<double> psigCB = {0.0136624, 0.230678, 1.73254};
        sigCB = sqrt(pow(psigCB[0],2) + pow(psigCB[1],2)/mll + pow(psigCB[2]/mll, 2));
        alpCB = 1.59112;
        static const vector<double> pnCB = {1.13055, 0.76705, 0.00298312};
        nCB = pnCB[0] + pnCB[1]*exp(-pnCB[2]*mll);
        static const vector<double> pmuG = {-0.00402708, 0.814172, -3.94281e-7, 7.97076e-6, -87.6397, -1.64806e-11};
        muG = pmuG[0] + pmuG[1]/mll + pmuG[2]*mll + pmuG[3]*pow(lnm,3) + pmuG[4]/sqr(mll) + pmuG[5]*sqr(mll);
        static const vector<double> psigG = {0.00553858, 0.140909, 0.644418};
        sigG = sqrt(sqr(psigG[0]) + sqr(psigG[1])/mll + sqr(psigG[2]/mll));
        static const vector<double> pkappa = {0.347003, 0.135768, 0.00372497, -2.2822e-5, 5.06351e-13};
        kappa = pkappa[0] + pkappa[1]*exp(-pkappa[2]*mll) + pkappa[3]*mll + pkappa[4]*pow(mll,3);
      } else {
        static const vector<double> pmuCB = {-0.0891397, 10.6169, -951.712, 74775.3, 5.60192e-5, -1.58827e-9, -3.81706e-13};
        muCB = pmuCB[0] + pmuCB[1]/mll + pmuCB[2]/sqr(mll) + pmuCB[3]/pow(mll,3) + pmuCB[4]*mll + pmuCB[5]*sqr(mll) + pmuCB[6]*pow(mll,3);
        static const vector<double> psigCB = {0.0836349, -8.98476, 491.19, 5.18068e-5, -3.45042e-10};
        sigCB = psigCB[0] + psigCB[1]/mll + psigCB[2]/sqr(mll) + psigCB[3]*mll + psigCB[4]*sqr(mll);
        static const vector<double> palpCB = {0.512577, 252.922, -79337.4, 7.31863e6, 0.000237883};
        alpCB = palpCB[0] + palpCB[1]/mll + palpCB[2]/sqr(mll) + palpCB[3]/pow(mll,3) + palpCB[4]*mll;
        static const vector<double> pnCB = {1.13055, 0.76705, 0.00298312};
        nCB = 6.08818;
        static const vector<double> pmuG = {-0.00410659, 2.82352e-6, -6.264e-9, 1.25547e-12, -9.94431e-17};
        muG = pmuG[0] + pmuG[1]*mll + pmuG[2]*sqr(mll) + pmuG[3]*pow(mll,3) + pmuG[4]*pow(mll,4);
        static const vector<double> psigG = {0.0214264, -0.795058, 15.4726, 3.38205e-5, -1.64984e-9};
        sigG = psigG[0] + psigG[1]/mll + psigG[2]/sqr(mll) + psigG[3]*mll + psigG[4]*sqr(mll);
        static const vector<double> pkappa = {0.235477, -31.2227, 3447.34, 4.54408e-5, -3.25374e-9};
        kappa = pkappa[0] + pkappa[1]/mll + pkappa[2]/sqr(mll) + pkappa[3]*mll + pkappa[4]*sqr(mll);
      }

      // Do the smearing
      double mll_scale = -1;
      do {
        mll_scale = (rand01() > kappa) ? randnorm(muG, sigG) : randcrystalball(alpCB, nCB, muCB, sigCB);
      } while (fabs(mll_scale) > 0.5);
      const double mll_reco = (1 + mll_scale) * mll;

      // Require high-Mll to avoid the Z peak
      if (mll_reco < 225*GeV) vetoEvent;

      // Fill Mll histograms
      (mumu ? _h_mmm : _h_mee)->fill(mll_reco/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_mee, crossSection()*luminosity()/femtobarn/sumOfWeights());
      scale(_h_mmm, crossSection()*luminosity()/femtobarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_mee, _h_mmm;
    //@}


  };


  DECLARE_RIVET_PLUGIN(ATLAS_2019_I1725190);

}

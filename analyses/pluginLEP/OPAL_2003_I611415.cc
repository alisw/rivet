// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/GammaGammaFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Dijet production in photon-photon collisions at 198 GeV
  class OPAL_2003_I611415 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_2003_I611415);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // get the hadronic final state
      const GammaGammaKinematics& diskin = declare(GammaGammaKinematics(), "Kinematics");
      const FinalState & fs = declare(GammaGammaFinalState(diskin), "FS");
      declare(FastJets(fs, FastJets::KT,1.),"Jets");
      book(_h_theta[0]   ,  1,1,1);
      book(_h_theta[1]   ,  2,1,1);
      book(_h_ET[0]      ,  3,1,1);
      book(_h_ET[1]      ,  4,1,1);
      book(_h_ET[2]      ,  5,1,1);
      book(_h_xg[0][0]   ,  6,1,1);
      book(_h_xg[0][1]   ,  7,1,1);
      book(_h_xg[1][0]   ,  9,1,1);
      book(_h_xg[1][1]   , 10,1,1);
      book(_h_xg[2][0]   , 11,1,1);
      book(_h_xg[2][1]   , 12,1,1);
      book(_h_xg_high    ,  8,1,1);
      book(_h_xlog[0]    , 13,1,1);
      book(_h_xlog[1]    , 14,1,1);
      book(_h_xlog[2]    , 15,1,1);
      book(_h_eta_diff[0], 16,1,1);
      book(_h_eta_diff[1], 17,1,1);
      book(_h_eta_min[0] , 18,1,1);
      book(_h_eta_min[1] , 19,1,1);
      book(_h_eta_max[0] , 20,1,1);
      book(_h_eta_max[1] , 21,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // need at least two jets with |eta|<2 and pT>3
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::Et > 3.*GeV and Cuts::abseta < 2.);
      if(jets.size()<2) vetoEvent;
      if(jets[0].Et()<jets[1].Et()) swap(jets[0],jets[1]);
      // Ets of jets
      double Et1 = jets[0].Et(), Et2 = jets[1].Et();
      // average Et
      double Etbar = 0.5*(Et1+Et2);
      double etaBar = 0.5*(jets[0].eta()+jets[1].eta());
      if(Etbar<5.) vetoEvent;
      // assymetry cut
      if((Et1-Et2)/(Et1+Et2)>.25) vetoEvent;
      // calculate x_gamma
      FourMomentum psum;
      for(const Particle & part : apply<FinalState>(event,"FS").particles()) {
        psum += part.momentum();
      }
      FourMomentum pj = jets[0].momentum()+jets[1].momentum();
      double xp = (pj.E()+pj.pz())/(psum.E()+psum.pz());
      double xm = (pj.E()-pj.pz())/(psum.E()-psum.pz());
      double cost = tanh(0.5*(jets[0].eta()-jets[1].eta()));
      // cost distributions
      if(pj.mass()>15.*GeV && etaBar<=1.) {
        if(xp>0.75 && xm>0.75)
          _h_theta[0]->fill(abs(cost));
        else if(xp<0.75 && xm<0.75)
          _h_theta[1]->fill(abs(cost));
      }
      // ET distributions
      _h_ET[0]->fill(Etbar);
      if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75))
        _h_ET[1]->fill(Etbar);
      else if(xp<0.75 && xm <0.75)
        _h_ET[2]->fill(Etbar);
      if(Etbar>=5.&&Etbar<7.) {
        _h_xg[0][0]->fill(xp);
        _h_xg[0][0]->fill(xm);
        _h_xlog[0]->fill(log(xp));
        _h_xlog[0]->fill(log(xm));
        if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75)) {
          _h_xg[1][0]->fill(xp);
          _h_xg[1][0]->fill(xm);
          _h_xlog[1]->fill(log(xp));
          _h_xlog[1]->fill(log(xm));
        }
        else if(xp<0.75 && xm <0.75) {
          _h_xg[2][0]->fill(xp);
          _h_xg[2][0]->fill(xm);
          _h_xlog[2]->fill(log(xp));
          _h_xlog[2]->fill(log(xm));
        }
      }
      else if(Etbar>=7.&& Etbar<11.) {
        _h_xg[0][1]->fill(xp);
        _h_xg[0][1]->fill(xm);
        if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75)) {
          _h_xg[1][1]->fill(xp);
          _h_xg[1][1]->fill(xm);
        }
        else if(xp<0.75 && xm <0.75) {
          _h_xg[2][1]->fill(xp);
          _h_xg[2][1]->fill(xm);
        }
      }
      else if(Etbar>=11.&& Etbar<25.) {
        _h_xg_high->fill(xp);
        _h_xg_high->fill(xm);
      }
      // vs eta
      double etaMin = min(abs(jets[0].eta()),abs(jets[1].eta()));
      double etaMax = max(abs(jets[0].eta()),abs(jets[1].eta()));
      if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75)) {
        _h_eta_diff[0]->fill(abs(jets[0].eta()-jets[1].eta()));
        _h_eta_min[0]->fill(etaMin);
        _h_eta_max[0]->fill(etaMax);
      }
      else if(xp<0.75 && xm <0.75) {
        _h_eta_diff[1]->fill(abs(jets[0].eta()-jets[1].eta()));
        _h_eta_min[1]->fill(etaMin);
        _h_eta_max[1]->fill(etaMax);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_h_theta[ix], fact);
        scale(_h_eta_diff[ix], fact);
        scale(_h_eta_min[ix], fact);
        scale(_h_eta_max[ix], fact);
        for(unsigned int iy=0;iy<3;++iy) {
          scale(_h_xg[iy][ix],fact);
        }
      }
      for(unsigned int ix=0;ix<3;++ix) {
        scale(_h_ET[ix],fact);
        scale(_h_xlog[ix],fact);
      }
      scale(_h_xg_high,fact);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_theta[2],_h_ET[3],_h_xg[3][2],_h_xg_high;
    Histo1DPtr _h_xlog[3],_h_eta_diff[2],_h_eta_min[2],_h_eta_max[2];
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(OPAL_2003_I611415);


}

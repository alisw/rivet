// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Multi-jet cross-sections in charged current $e^{\pm} p$ scattering at HERA
  class ZEUS_2008_I780108 : public Analysis {
  public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_2008_I780108);

      /// @name Analysis methods
      /// @{

      /// Book histograms and initialise projections before the run
      void init() {

          // Projections
          const DISKinematics diskin;
          declare(diskin, "Kinematics");
          const DISFinalState disfs(DISFinalState::BoostFrame::LAB);
          FastJets jets(disfs, FastJets::KT, 1.0);
          declare(jets, "Jets");

          // Table 11
          book(_h_eta_incl[0], 11, 1, 1);
          book(_h_eta_incl[1], 11, 1, 2);
          // Table 12
          book(_h_eta_di[0], 12, 1, 1);
          book(_h_eta_di[1], 12, 1, 2);
          // Table 13
          book(_h_eta_tri[0], 13, 1, 1);
          book(_h_eta_tri[1], 13, 1, 2);
          // Table 14
          book(_h_et_incl[0], 14, 1, 1);
          book(_h_et_incl[1], 14, 1, 2);
          // Table 15
          book(_h_et_di[0], 15, 1, 1);
          book(_h_et_di[1], 15, 1, 2);
          // Table 16
          book(_h_et_tri[0], 16, 1, 1);
          book(_h_et_tri[1], 16, 1, 2);
          // Table 17
          book(_h_q2_incl[0], 17, 1, 1);
          book(_h_q2_incl[1], 17, 1, 2);
          // Table 18
          book(_h_q2_di[0], 18, 1, 1);
          book(_h_q2_di[1], 18, 1, 2);
          // Table 19
          book(_h_q2_tri[0], 19, 1, 1);
          book(_h_q2_tri[1], 19, 1, 2);

          // Table 20
          book(_h_x_incl[0], 20, 1, 1);
          book(_h_x_incl[1], 20, 1, 2);
          // Table 22
          book(_h_m_di[0], 22, 1, 1);
          book(_h_m_di[1], 22, 1, 2);
          // Table 23
          book(_h_m_tri[0], 23, 1, 1);
          book(_h_m_tri[1], 23, 1, 2);

      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

          int fLepton = 0;

          const ParticlePair bs = event.beams();
          if (bs.first.pid() == PID::POSITRON || bs.second.pid() == PID::POSITRON) fLepton = 1;
          const Particle& bproton = (bs.first.pid() == PID::PROTON) ? bs.first : bs.second;
          const int orientation = sign(bproton.momentum().pz());
          // DIS kinematics
          const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
          double q2  = dk.Q2();
          double x   = dk.x();
          double y   = dk.y();

          if (q2 < 200) vetoEvent;
          if (y > 0.9) vetoEvent;
          // Jet selection
          const Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 5*GeV && Cuts::etaIn(-1*orientation, 2.5*orientation), cmpMomByEt);
          MSG_DEBUG("Jet multiplicity = " << jets.size());
          if (jets.size() < 1) vetoEvent;
          if (jets[0].Et() < 14*GeV) vetoEvent;

          double eta12 = 0;
          double et12 = 0;

          double et123 = 0;
          double eta123 =0;

          for (size_t i = 0; i < jets.size(); i++)
          {
              if (jets[i].Et() < 14*GeV) continue;
              _h_eta_incl[fLepton]->fill(orientation*jets[i].eta());
              _h_et_incl[fLepton]->fill(jets[i].Et());
              _h_q2_incl[fLepton]->fill(q2);
              _h_x_incl[fLepton]->fill(x);
          }

          if (jets.size() > 1)
          {
              eta12 = orientation*(jets[0].eta() + jets[1].eta())/2;
              et12 = (jets[0].Et() + jets[1].Et())/2;
              _h_eta_di[fLepton]->fill(eta12);
              _h_et_di[fLepton]->fill(et12);
              _h_q2_di[fLepton]->fill(q2);
              _h_m_di[fLepton]->fill(   (jets[0].momentum()+jets[1].momentum()).mass());
          }

          if (jets.size() > 2)
          {
              eta123 = orientation*(jets[0].eta() + jets[1].eta()+jets[2].eta())/3;
              et123 = (jets[0].Et() + jets[1].Et()+jets[2].Et())/3;

              _h_eta_tri[fLepton]->fill(eta123);
              _h_et_tri[fLepton]->fill(et123);
              _h_q2_tri[fLepton]->fill(q2);
              _h_m_tri[fLepton]->fill(   (jets[0].momentum()+jets[1].momentum()+jets[2].momentum()).mass());
          }
      }


      /// Normalise histograms etc., after the run
      void finalize() {
          const double sf = crossSection()/picobarn/sumOfWeights();


          scale(_h_eta_incl, sf);
          scale(_h_eta_di,   sf);
          scale(_h_eta_tri,  sf);
          scale(_h_et_incl,  sf);
          scale(_h_et_di,    sf);
          scale(_h_et_tri,   sf);
          scale(_h_q2_incl,  sf);
          scale(_h_q2_di,    sf);
          scale(_h_q2_tri,   sf);
          scale(_h_x_incl,   sf);
          scale(_h_m_di,     sf);
          scale(_h_m_tri,    sf);

      }

      /// @}


      /// @name Histograms
      /// @{
      Histo1DPtr _h_eta_incl[2], _h_eta_di[2], _h_eta_tri[2];
      Histo1DPtr _h_et_incl[2], _h_et_di[2], _h_et_tri[2];
      Histo1DPtr _h_q2_incl[2], _h_q2_di[2], _h_q2_tri[2];
      Histo1DPtr _h_x_incl[2], _h_x_di[2], _h_x_tri[2];
      Histo1DPtr _h_m_di[2], _h_m_tri[2];
      /// @}
  };


  RIVET_DECLARE_PLUGIN(ZEUS_2008_I780108);

}

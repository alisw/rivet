// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"


/// @author Peter Wijeratne <paw@hep.ucl.ac.uk>
/// @author Robindra Prabhu <prabhu@cern.ch>
namespace Rivet {

  // A very basic analysis sensitive to ET flow in minbias and dijet events
  class ATLAS_2012_I1183818 : public Analysis {

  public:

    ATLAS_2012_I1183818()
      : Analysis("ATLAS_2012_I1183818")
    {}


  public:

    void init() {

      const FinalState cnfs((Cuts::etaIn(-4.8, 4.8)));
      const ChargedFinalState cfs((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  250*MeV));
      declare(cnfs, "FS");
      declare(cfs, "CFS");

      const FastJets jetsAntiKt4(cnfs, FastJets::ANTIKT, 0.4);
      declare(jetsAntiKt4, "AntiKt4Jets");

      // ------- MINBIAS HISTOGRAMS --------
      //
      // MB event counter
      book(m_chargedEvents, "m_chargedEvents");

      book(_h_ETflowEta ,1, 1, 1);
      book(_h_SumETbin1 ,3, 1, 1);
      book(_h_SumETbin2 ,4, 1, 1);
      book(_h_SumETbin3 ,5, 1, 1);
      book(_h_SumETbin4 ,6, 1, 1);
      book(_h_SumETbin5 ,7, 1, 1);
      book(_h_SumETbin6 ,8, 1, 1);

      // ------- DIJET HISTOGRAMS --------
      //
      // Dijet event counter
      book(m_events_dijets, "m_chargedEvents");

      // sumET
      book(_h_transETflowEta , 2, 1, 1);
      book(_h_transSumETbin1 , 9, 1, 1);
      book(_h_transSumETbin2 ,10, 1, 1);
      book(_h_transSumETbin3 ,11, 1, 1);
      book(_h_transSumETbin4 ,12, 1, 1);
      book(_h_transSumETbin5 ,13, 1, 1);
      book(_h_transSumETbin6 ,14, 1, 1);


    }


    void analyze(const Event& event) {

      const FinalState& cfs = apply<FinalState>(event, "CFS");

      bool isCharged = false;
      if (cfs.size() >= 2) {  // event selection: > 2 charged particles with pT > 250.MeV and |eta| < 2.5
        isCharged = true;
        m_chargedEvents->fill();
      }

      const FinalState& cnfs = apply<FinalState>(event, "FS");

      Particles particles;
      for( const Particle& p : cnfs.particles() ) {
        // enforce truth selection representing detected particle sensitivity
        double pp = p.p3().mod();
        if (PID::charge3(p.pid()) != 0 && pp < 0.5*GeV) continue;
        if (PID::charge3(p.pid()) == 0 && pp < 0.2*GeV) continue;

        particles.push_back(p);
      }


      // get jets
      const FastJets& jetsAntiKt4 = apply<FastJets>(event, "AntiKt4Jets");
      const Jets& jets = jetsAntiKt4.jetsByPt(20.0*GeV);

      // initialise sumET variables
      double sumETbin1 = 0;
      double sumETbin2 = 0;
      double sumETbin3 = 0;
      double sumETbin4 = 0;
      double sumETbin5 = 0;
      double sumETbin6 = 0;

      // if (passes event selection)
      if (isCharged) {

        for( const Particle& p : particles ) {

          ///calculate variables
          double ET = p.Et()/GeV;
          double eta = p.abseta();

          // fill histograms
          _h_ETflowEta->fill(eta, ET);

          if      (eta <  0.8) sumETbin1 += ET;
          else if (eta <  1.6) sumETbin2 += ET;
          else if (eta <  2.4) sumETbin3 += ET;
          else if (eta <  3.2) sumETbin4 += ET;
          else if (eta <  4.0) sumETbin5 += ET;
          else if (eta <= 4.8) sumETbin6 += ET;

        } // end of for

        _h_SumETbin1->fill(sumETbin1);
        _h_SumETbin2->fill(sumETbin2);
        _h_SumETbin3->fill(sumETbin3);
        _h_SumETbin4->fill(sumETbin4);
        _h_SumETbin5->fill(sumETbin5);
        _h_SumETbin6->fill(sumETbin6);
      }

      // --- do dijet analysis ---

      if ( jets.size() >= 2                       && // require at least two jets
           jets[0].Et() >= 20.*GeV     && // require two leading jets to pass ET cuts
           jets[1].Et() >= 20.*GeV     &&
           fabs(jets[0].eta()) < 2.5   && // require leading jets to be central
           fabs(jets[1].eta()) < 2.5   &&
           deltaPhi(jets[0], jets[1]) > 2.5       && // require back-to-back topology
           jets[1].Et()/jets[0].Et() >= 0.5) { //require ET-balance

        // found an event that satisfies dijet selection, now fill histograms...
        // initialise dijet sumET variables
        double trans_sumET_bin1 = 0.;
        double trans_sumET_bin2 = 0.;
        double trans_sumET_bin3 = 0.;
        double trans_sumET_bin4 = 0.;
        double trans_sumET_bin5 = 0.;
        double trans_sumET_bin6 = 0.;

        m_events_dijets->fill();

        // loop over all particles and check their relation to leading jet
        for( const Particle& particle : particles ) {

          // calculate variables
          double dPhi = deltaPhi( jets[0], particle.momentum() );
          double ET   = particle.Et()/GeV;
          double eta  = fabs(particle.eta());

          // Transverse region
          if ( dPhi > 1./3.*M_PI && dPhi < 2./3.*M_PI ) {
            _h_transETflowEta->fill( eta, ET );
            if      (eta <  0.8) { trans_sumET_bin1 += ET; }
            else if (eta <  1.6) { trans_sumET_bin2 += ET; }
            else if (eta <  2.4) { trans_sumET_bin3 += ET; }
            else if (eta <  3.2) { trans_sumET_bin4 += ET; }
            else if (eta <  4.0) { trans_sumET_bin5 += ET; }
            else if (eta <= 4.8) { trans_sumET_bin6 += ET; }
          }

        } // end loop over particles

        _h_transSumETbin1->fill( trans_sumET_bin1);
        _h_transSumETbin2->fill( trans_sumET_bin2);
        _h_transSumETbin3->fill( trans_sumET_bin3);
        _h_transSumETbin4->fill( trans_sumET_bin4);
        _h_transSumETbin5->fill( trans_sumET_bin5);
        _h_transSumETbin6->fill( trans_sumET_bin6);
      } // end of dijet selection cuts

    }


    void finalize() {
      /// several scale factors here:
      /// 1. nEvents (m_chargedEvents)
      /// 2. phase-space (2*M_PI)
      /// 3. double binning due to symmetrisation (2)
      scale( _h_ETflowEta, 1./m_chargedEvents->val()/(4.*M_PI) );
      scale( _h_SumETbin1, 1./m_chargedEvents->val() );
      scale( _h_SumETbin2, 1./m_chargedEvents->val() );
      scale( _h_SumETbin3, 1./m_chargedEvents->val() );
      scale( _h_SumETbin4, 1./m_chargedEvents->val() );
      scale( _h_SumETbin5, 1./m_chargedEvents->val() );
      scale( _h_SumETbin6, 1./m_chargedEvents->val() );

      //Dijet analysis

      // Dijet scale factors:
      //1. number of events passing dijet selection
      //2. phase-space: 1. / 2/3*M_PI
      //3. double binning due to symmetrisation in |eta| plots : 1/2
      scale( _h_transETflowEta, 1./m_events_dijets->val() * 1./(4./3.*M_PI) );
      scale( _h_transSumETbin1, 1./m_events_dijets->val() );
      scale( _h_transSumETbin2, 1./m_events_dijets->val() );
      scale( _h_transSumETbin3, 1./m_events_dijets->val() );
      scale( _h_transSumETbin4, 1./m_events_dijets->val() );
      scale( _h_transSumETbin5, 1./m_events_dijets->val() );
      scale( _h_transSumETbin6, 1./m_events_dijets->val() );
    }

  private:

    // Event counts
    CounterPtr m_chargedEvents;
    CounterPtr m_events_dijets;

    // Minbias-analysis: variable + histograms
    Histo1DPtr _h_ETflowEta;
    Histo1DPtr _h_SumETbin1;
    Histo1DPtr _h_SumETbin2;
    Histo1DPtr _h_SumETbin3;
    Histo1DPtr _h_SumETbin4;
    Histo1DPtr _h_SumETbin5;
    Histo1DPtr _h_SumETbin6;

    // Transverse region
    Histo1DPtr _h_transETflowEta;
    Histo1DPtr _h_transSumETbin1;
    Histo1DPtr _h_transSumETbin2;
    Histo1DPtr _h_transSumETbin3;
    Histo1DPtr _h_transSumETbin4;
    Histo1DPtr _h_transSumETbin5;
    Histo1DPtr _h_transSumETbin6;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1183818);

}

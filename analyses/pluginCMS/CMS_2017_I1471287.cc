// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Tools/Correlators.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2017_I1471287 : public CumulantAnalysis {
  public:

    /// Constructor
    CMS_2017_I1471287() : CumulantAnalysis("CMS_2017_I1471287") {
    
    };


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // A projection for charged tracks to manage centrality, corresponding
      // to CMS offline tracks.
      ChargedFinalState cfsMult(Cuts::abseta < 2.4 && Cuts::pT > 0.4*GeV);
      addProjection(cfsMult, "CFSMult");
      
      // The positive eta side used for rapidity gap, integrated.
      const ChargedFinalState& cfsp = ChargedFinalState(Cuts::eta > 1.0 && 
        Cuts::eta < 2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 3.0*GeV);
      declare(cfsp, "CFSP");
      // ..negative ditto.
      const ChargedFinalState& cfsn = ChargedFinalState(Cuts::eta < -1.0 && 
        Cuts::eta > -2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 3.0*GeV);
      declare(cfsn, "CFSN");


      // The positive eta side used for rapidity gap, differential, charged particles.
      const ChargedFinalState& cfsppT = ChargedFinalState(Cuts::eta > 1.0 && 
        Cuts::eta < 2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 6.0*GeV);
      declare(cfsppT, "CFSPPT");
      // ..negative ditto.
      const ChargedFinalState& cfsnpT = ChargedFinalState(Cuts::eta < -1.0 && 
        Cuts::eta > -2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 6.0*GeV);
      declare(cfsnpT, "CFSNPT");

      // The positive eta side used for rapidity gap, differential, Kaons.
      const PrimaryParticles& kfsppT = PrimaryParticles({310},Cuts::eta > 1.0 && 
        Cuts::eta < 2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 6.0*GeV);
      declare(kfsppT, "KFSP");
      // ..negative ditto.
      const PrimaryParticles& kfsnpT = PrimaryParticles({310},Cuts::eta < -1.0 && 
        Cuts::eta > -2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 6.0*GeV);
      declare(kfsnpT, "KFSN");
     // The positive eta side used for rapidity gap, differential, Lambda.
      const PrimaryParticles& lfsppT = PrimaryParticles({3122},Cuts::eta > 1.0 && 
        Cuts::eta < 2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 6.0*GeV);
      declare(lfsppT, "LFSP");
      // ..negative ditto.
      const PrimaryParticles& lfsnpT = PrimaryParticles({3122},Cuts::eta < -1.0 && 
        Cuts::eta > -2.0 && Cuts::pT > 0.3*GeV && Cuts::pT < 6.0*GeV);
      declare(lfsnpT, "LFSN");
      
      // v22 |delta eta| > 2 (fig 4a)
      h_v22 = bookScatter2D(1, 1, 1, true);
      // v32 |delta eta| > 2 (fig 4b)
      h_v32 = bookScatter2D(3, 1, 1, true);
      // v22(pT) high mult., high pT (fig 6a)
      h_v22pT = bookScatter2D(11, 1, 1, true);
      // v22(pT) charged low mult. (fig. 7a)
      h_v22pTh = bookScatter2D(17, 1, 1, true);
      // v22(pT) K0S low mult. (fig. 7a)
      h_v22pTK = bookScatter2D(18, 1, 1, true);
      // v22(pT) Lambda low mult. (fig. 7a)
      h_v22pTL = bookScatter2D(19, 1, 1, true);
      // v22(pT) K0S high mult. (fig. 7b)
      h_v22pTKc = bookScatter2D(21, 1, 1, true);
      // v22(pT) Lambda high mult. (fig. 7b)
      h_v22pTLc = bookScatter2D(22, 1, 1, true);
      // c24 (fig. 9a)
      h_c24 = bookScatter2D(28, 1, 1, true);
      // c26 (fig. 9b)
      h_c26 = bookScatter2D(31, 1, 1, true);

      // Corresponding event averaged correlators.
      ec22 = bookECorrelatorGap<2,2>("ec22",h_v22);
      ec32 = bookECorrelatorGap<3,2>("ec32",h_v32);

      // ... pT binned
      ec22pT = bookECorrelatorGap<2,2>("ec22pT",h_v22pT);
      ec22pTh = bookECorrelatorGap<2,2>("ec22pTh",h_v22pTh);
      ec22pTK = bookECorrelatorGap<2,2>("ec22pTK",h_v22pTK);
      ec22pTL = bookECorrelatorGap<2,2>("ec22pTL",h_v22pTL);
      ec22pTKc = bookECorrelatorGap<2,2>("ec22pTKc",h_v22pTKc);
      ec22pTLc = bookECorrelatorGap<2,2>("ec22pTLc",h_v22pTLc);

      // Maximal N and P for the gapped.
      pair<int, int> max = getMaxValues(); 
      
      // For the four particle cumulant.
      ec22_4 = bookECorrelator<2,2>("ec22_4",h_c24);
      ec24_4 = bookECorrelator<2,4>("ec24_4",h_c24);

      // For the six particle cumulant.
      ec22_6 = bookECorrelator<2,2>("ec22_6",h_c26);
      ec24_6 = bookECorrelator<2,4>("ec24_6",h_c26);
      ec26_6 = bookECorrelator<2,6>("ec26_6",h_c26);

      // Maximal N and P for the higher orders.
      pair<int, int> maxH = getMaxValues();

      // Declare correlator projections.
      // For integrated.
      declare(Correlators(cfsMult, maxH.first, maxH.second),"CH");

      // ... gapped
      declare(Correlators(cfsp, max.first, max.second),"CPos");
      declare(Correlators(cfsn, max.first, max.second),"CNeg");

      // For pT differential, charged particles, low multiplicity.
      declare(Correlators(cfsppT, max.first, max.second, h_v22pTh),"CPosLowPT");
      declare(Correlators(cfsnpT, max.first, max.second, h_v22pTh),"CNegLowPT");

      // For pT differential, charged particles, high multiplicity.
      declare(Correlators(cfsppT, max.first, max.second, h_v22pT),"CPosHighPT");
      declare(Correlators(cfsnpT, max.first, max.second, h_v22pT),"CNegHighPT");
      
      // For pT differential, kaons. low multiplicity.
      declare(Correlators(kfsppT, max.first, max.second, h_v22pTK),"CPosLowPTK");
      declare(Correlators(kfsnpT, max.first, max.second, h_v22pTK),"CNegLowPTK");
     
      // For pT differential, kaons. high multiplicity.
      declare(Correlators(kfsppT, max.first, max.second, h_v22pTKc),"CPosHighPTK");
      declare(Correlators(kfsnpT, max.first, max.second, h_v22pTKc),"CNegHighPTK");
      
      // For pT differential, lambda. low multiplicity.
      declare(Correlators(lfsppT, max.first, max.second, h_v22pTL),"CPosLowPTL");
      declare(Correlators(lfsnpT, max.first, max.second, h_v22pTL),"CNegLowPTL");
      
      // For pT differential, lambda. high multiplicity.
      declare(Correlators(lfsppT, max.first, max.second, h_v22pTLc),"CPosHighPTL");
      declare(Correlators(lfsnpT, max.first, max.second, h_v22pTLc),"CNegHighPTL");
      

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double w = event.weight();
      const double nTrk = apply<ChargedFinalState>(event, "CFSMult").particles().size();
      if (nTrk < 10) vetoEvent;

      // The correlators.
      const Correlators& ch = apply<Correlators>(event, "CH");

      const Correlators& cp = apply<Correlators>(event, "CPos");
      const Correlators& cn = apply<Correlators>(event, "CNeg");

      const Correlators& cpLow = apply<Correlators>(event, "CPosLowPT");
      const Correlators& cnLow = apply<Correlators>(event, "CNegLowPT");

      const Correlators& cpHigh = apply<Correlators>(event, "CPosHighPT");
      const Correlators& cnHigh = apply<Correlators>(event, "CNegHighPT");

      const Correlators& cpLowK = apply<Correlators>(event, "CPosLowPTK");
      const Correlators& cnLowK = apply<Correlators>(event, "CNegLowPTK");

      const Correlators& cpHighK = apply<Correlators>(event, "CPosHighPTK");
      const Correlators& cnHighK = apply<Correlators>(event, "CNegHighPTK");
      
      const Correlators& cpLowL = apply<Correlators>(event, "CPosLowPTL");
      const Correlators& cnLowL = apply<Correlators>(event, "CNegLowPTL");

      const Correlators& cpHighL = apply<Correlators>(event, "CPosHighPTL");
      const Correlators& cnHighL = apply<Correlators>(event, "CNegHighPTL");

      ec22->fill(nTrk, cp, cn, w);
      ec32->fill(nTrk, cp, cn, w);

      ec22_4->fill(nTrk, ch, w);
      ec24_4->fill(nTrk, ch, w);
      ec22_6->fill(nTrk, ch, w);
      ec24_6->fill(nTrk, ch, w);
      ec26_6->fill(nTrk, ch, w);

      if (nTrk < 20) {
        ec22pTh->fill(cpLow, cnLow, w);
        ec22pTK->fill(cpLowK, cnLowK, w);
        ec22pTL->fill(cpLowL, cnLowL, w);
      }
      else if(nTrk >= 105 && nTrk < 150)
        ec22pT->fill(cpHigh, cnHigh, w);
        ec22pTKc->fill(cpHighK, cnHighK, w);
        ec22pTLc->fill(cpHighL, cnHighL, w);
      }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Correlators must be streamed 
      // in order to run reentrant finalize.
      stream();
      cnTwoInt(h_v22, ec22);
      cnTwoInt(h_v32, ec32);
      vnTwoDiff(h_v22pT, ec22pT);
      vnTwoDiff(h_v22pTh, ec22pTh);
      cnFourInt(h_c24, ec22_4, ec24_4);
      cnSixInt(h_c26, ec22_6, ec24_6, ec26_6);

      // Set correct reference flow for pid flow.
      ec22pTK->setReference(ec22pTh->getReference());
      vnTwoDiff(h_v22pTK, ec22pTK);
      ec22pTL->setReference(ec22pTh->getReference());
      vnTwoDiff(h_v22pTL, ec22pTL);
      ec22pTKc->setReference(ec22pT->getReference());
      vnTwoDiff(h_v22pTKc, ec22pTKc);
      ec22pTLc->setReference(ec22pT->getReference());
      vnTwoDiff(h_v22pTLc, ec22pTLc);

    }

    //@}
      Scatter2DPtr h_v22;
      Scatter2DPtr h_v32;
      Scatter2DPtr h_v22pT;
      Scatter2DPtr h_v22pTh;
      Scatter2DPtr h_v22pTK;
      Scatter2DPtr h_v22pTL;
      Scatter2DPtr h_v22pTKc;
      Scatter2DPtr h_v22pTLc;
      Scatter2DPtr h_c24;
      Scatter2DPtr h_c26;

      ECorrPtr ec22;
      ECorrPtr ec32;

      ECorrPtr ec22_4;
      ECorrPtr ec24_4;

      ECorrPtr ec22_6;
      ECorrPtr ec24_6;
      ECorrPtr ec26_6;

      ECorrPtr ec22pT;
      ECorrPtr ec22pTh;
      ECorrPtr ec22pTK;
      ECorrPtr ec22pTL;
      ECorrPtr ec22pTKc;
      ECorrPtr ec22pTLc;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1471287);


}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Correlators.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {


  /// @brief Multiparticle azimuthal correlations PbPb 5TeV.
  class ALICE_2016_I1419244 : public CumulantAnalysis {
  public:

    /// Constructor
    ALICE_2016_I1419244()
      : CumulantAnalysis("ALICE_2016_I1419244") {
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Declare the trigger projection.
      declare<ALICE::V0AndTrigger>(ALICE::V0AndTrigger(),"V0-AND");
      // Centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M","V0M");

      // The full central charged final state.
      const ChargedFinalState& cfs = ChargedFinalState(Cuts::abseta < 0.8 &&
        Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfs, "CFS");

      // The positive eta side used for rapidity gap.
      const ChargedFinalState& cfsp = ChargedFinalState(Cuts::eta > 0.5 &&
        Cuts::eta < 0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfsp, "CFSP");
      // ..negative ditto.
      const ChargedFinalState& cfsn = ChargedFinalState(Cuts::eta < -0.5 &&
        Cuts::eta > -0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfsn, "CFSN");

      // Book histograms before booking the correlators
      // to have access to bin edges.
      book(h_v22gap, 1, 1, 1, true);
      book(h_v24, 1, 1, 2, true);
      book(h_v26, 1, 1, 3, true);
      book(h_v28, 1, 1, 4, true);
      book(h_v32gap, 2, 1, 1, true);
      book(h_v42gap, 2, 1, 2, true);
      book(h_v22gappT, 8, 1, 1, true);
      book(h_v32gappT, 8, 1, 2, true);
      book(h_v42gappT, 8, 1, 3, true);
      book(h_v24pT10, 9, 1, 1, true);
      book(h_v24pT20, 9, 1, 2, true);
      book(h_v24pT30, 9, 1, 3, true);

      book(h_c22gap, "_c22gap", refData(1, 1, 1));
      book(h_c24, "_c24", refData(1, 1, 2));
      book(h_c26, "_c26", refData(1, 1, 3));
      book(h_c28, "_c28", refData(1, 1, 4));
      book(h_c32gap, "_c32gap", refData(8, 1, 2));
      book(h_c42gap, "_c24gap", refData(8, 1, 3));

      book(h_ec22gap, "_ec22gap", refData(1, 1, 1));
      book(h_ec22, "_ec22", refData(1, 1, 2));
      book(h_ec24, "_ec24", refData(1, 1, 2));
      book(h_ec26, "_ec26", refData(1, 1, 3));
      book(h_ec28, "_ec28", refData(1, 1, 4));

      // Corresponding event averaged correlators.
      // Integrated, with gap.
      ec22gap = bookECorrelatorGap<2,2>("ec22gap",refData(1,1,1));
      ec32gap = bookECorrelatorGap<3,2>("ec32gap",refData(2,1,1));
      ec42gap = bookECorrelatorGap<4,2>("ec42gap",refData(2,1,2));

      // Integrated, no gap.
      ec22 = bookECorrelator<2,2>("ec22",refData(1,1,2));
      ec24 = bookECorrelator<2,4>("ec24",refData(1,1,2));
      ec26 = bookECorrelator<2,6>("ec26",refData(1,1,3));
      ec28 = bookECorrelator<2,8>("ec28",refData(1,1,4));

      // pT differential, no gap, three centralities.
      ec22pT10 = bookECorrelator<2,2>("ec22pT10",refData(9,1,1));
      ec24pT10 = bookECorrelator<2,4>("ec24pT10",refData(9,1,1));

      ec22pT20 = bookECorrelator<2,2>("ec22pT20",refData(9,1,2));
      ec24pT20 = bookECorrelator<2,4>("ec24pT20",refData(9,1,2));

      ec22pT30 = bookECorrelator<2,2>("ec22pT30",refData(9,1,3));
      ec24pT30 = bookECorrelator<2,4>("ec24pT30",refData(9,1,3));

      // pT differential, with gap, 30-40% centrality.
      ec22gappT = bookECorrelatorGap<2,2>("ec22gappT",refData(8,1,1));
      ec32gappT = bookECorrelatorGap<3,2>("ec32gappT",refData(8,1,2));
      ec42gappT = bookECorrelatorGap<4,2>("ec42gappT",refData(8,1,3));

      pair<int, int> max = getMaxValues();
      // Declare correlator projections.
      declare(Correlators(cfs, max.first, max.second, refData(8,1,1)),
        "Correlators");
      declare(Correlators(cfsp, max.first, max.second, refData(8,1,1)),
        "CorrelatorsPos");
      declare(Correlators(cfsn, max.first, max.second, refData(8,1,1)),
        "CorrelatorsNeg");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;

      // The centrality projection.
      const CentralityProjection& centProj =
        apply<CentralityProjection>(event,"V0M");

      // The centrality.
      const double cent = centProj();

      // The correlators projections.
      const Correlators& c = applyProjection<Correlators>(event,"Correlators");
      const Correlators& cp =
        applyProjection<Correlators>(event,"CorrelatorsPos");
      const Correlators& cn =
        applyProjection<Correlators>(event,"CorrelatorsNeg");

      ec22gap->fill(cent, cp, cn);
      ec32gap->fill(cent, cp, cn);
      ec42gap->fill(cent, cp, cn);

      ec22->fill(cent, c);
      ec24->fill(cent, c);
      ec26->fill(cent, c);
      ec28->fill(cent, c);

      if (cent < 10.) ec22pT10->fill(c), ec24pT10->fill(c);
      else if (cent < 20.) ec22pT20->fill(c), ec24pT20->fill(c);
      else if (cent < 30.) ec22pT30->fill(c), ec24pT30->fill(c);
      else if (cent < 40.) {
        ec22gappT->fill(cp, cn);
        ec32gappT->fill(cp, cn);
        ec42gappT->fill(cp, cn);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
       // Filling test histos.
       cnTwoInt(h_c22gap, ec22gap);
       cnTwoInt(h_c32gap, ec32gap);
       cnTwoInt(h_c42gap, ec42gap);
       cnFourInt(h_c24, ec22, ec24);
       cnSixInt(h_c26, ec22, ec24, ec26);
       cnEightInt(h_c28, ec22, ec24, ec26, ec28);

       corrPlot(h_ec22gap, ec22gap);
       corrPlot(h_ec22, ec22);
       corrPlot(h_ec24, ec24);
       corrPlot(h_ec26, ec26);
       corrPlot(h_ec28, ec28);

       vnTwoInt(h_v22gap, ec22gap);
       vnTwoInt(h_v32gap, ec32gap);
       vnTwoInt(h_v42gap, ec42gap);

       vnFourInt(h_v24, ec22, ec24);
       vnSixInt(h_v26, ec22, ec24, ec26);
       vnEightInt(h_v28, ec22, ec24, ec26, ec28);

       vnTwoDiff(h_v22gappT, ec22gappT);
       vnTwoDiff(h_v32gappT, ec32gappT);
       vnTwoDiff(h_v42gappT, ec42gappT);

       vnFourDiff(h_v24pT10, ec22pT10, ec24pT10);
       vnFourDiff(h_v24pT20, ec22pT20, ec24pT20);
       vnFourDiff(h_v24pT30, ec22pT30, ec24pT30);
       }


    //@}


    /// @name Histograms
    //@{
    // The integrated centrality dependent v2{n}.
    Scatter2DPtr h_v22gap;
    Scatter2DPtr h_v24;
    Scatter2DPtr h_v26;
    Scatter2DPtr h_v28;
    // The integrated, centrality dependent v3,4{2} gapped.
    Scatter2DPtr h_v32gap;
    Scatter2DPtr h_v42gap;
    // The pT differential, v2{2} gapped for 30-40% centrality
    Scatter2DPtr h_v22gappT;
    // ...v3{2} ditto.
    Scatter2DPtr h_v32gappT;
    // ... v4{2} ditto.
    Scatter2DPtr h_v42gappT;
    // The pT differential, centrality dependent v2{4}
    Scatter2DPtr h_v24pT10;
    Scatter2DPtr h_v24pT20;
    Scatter2DPtr h_v24pT30;

    // Test histograms -- cumulants
    Scatter2DPtr h_c22gap;
    Scatter2DPtr h_c24;
    Scatter2DPtr h_c26;
    Scatter2DPtr h_c28;
    Scatter2DPtr h_c32gap;
    Scatter2DPtr h_c42gap;

    // Test histograms -- correlators.
    Scatter2DPtr h_ec22gap;
    Scatter2DPtr h_ec22;
    Scatter2DPtr h_ec24;
    Scatter2DPtr h_ec26;
    Scatter2DPtr h_ec28;


    // The all event averaged correlators.
    // Integrated with gap.
    ECorrPtr ec22gap;
    ECorrPtr ec32gap;
    ECorrPtr ec42gap;

    // Integrated without gap.
    ECorrPtr ec22;
    ECorrPtr ec24;
    ECorrPtr ec26;
    ECorrPtr ec28;

    // pT dependent, 2 particle, gapped, 30-40% centrality
    ECorrPtr ec22gappT;
    ECorrPtr ec32gappT;
    ECorrPtr ec42gappT;

    // pT dependent, 4 particle, three centralities.
    ECorrPtr ec22pT10;
    ECorrPtr ec24pT10;

    ECorrPtr ec22pT20;
    ECorrPtr ec24pT20;

    ECorrPtr ec22pT30;
    ECorrPtr ec24pT30;


    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2016_I1419244);


}

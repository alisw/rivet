// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Correlators.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "YODA/Scatter2D.h"

namespace Rivet {

  /// @brief Multiparticle azimuthal correlations PbPb 5TeV.
  class ALICE_2016_I1419244 : public CumulantAnalysis {
  public:

    /// Constructor
    ALICE_2016_I1419244() : CumulantAnalysis("ALICE_2016_I1419244"){
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      // Initialise and register projections
      // Declare the trigger projection.
      declare<ALICE::V0AndTrigger>(ALICE::V0AndTrigger(),"V0-AND");
      // Centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality",
        "V0M","V0M");
      
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
      h_v22gap = bookScatter2D(1, 1, 1, true);
      h_v24 = bookScatter2D(1, 1, 2, true);
      h_v26 = bookScatter2D(1, 1, 3, true);
      h_v28 = bookScatter2D(1, 1, 4, true);
      h_v32gap = bookScatter2D(2, 1, 1, true);
      h_v42gap = bookScatter2D(2, 1, 2, true);
      h_v22gappT = bookScatter2D(8, 1, 1, true);
      h_v32gappT = bookScatter2D(8, 1, 2, true);
      h_v42gappT = bookScatter2D(8, 1, 3, true);
      h_v24pT10 = bookScatter2D(9, 1, 1, true);
      h_v24pT20 = bookScatter2D(9, 1, 2, true);
      h_v24pT30 = bookScatter2D(9, 1, 3, true);
     
      h_c22gap = bookScatter2D(h_v22gap,"/" + name() + "/c22gap");
      h_c24 = bookScatter2D(h_v24,"/" + name() + "/c24");
      h_c26 = bookScatter2D(h_v26,"/" + name() + "/c26");
      h_c28 = bookScatter2D(h_v28,"/" + name() + "/c28");
      h_c32gap = bookScatter2D(h_v32gap,"/" + name() + "/c32gap");
      h_c42gap = bookScatter2D(h_v42gap,"/" + name() + "/c42gap");

      h_ec22gap = bookScatter2D(h_v22gap,"/" + name() + "/ec22gap");
      h_ec22 = bookScatter2D(h_v24,"/" + name() + "/ec22");
      h_ec24 = bookScatter2D(h_v24,"/" + name() + "/ec24");
      h_ec26 = bookScatter2D(h_v26,"/" + name() + "/ec26");
      h_ec28 = bookScatter2D(h_v28,"/" + name() + "/ec28");

      
      // Corresponding event averaged correlators.
      // Integrated, with gap.
      ec22gap = bookECorrelatorGap<2,2>("ec22gap",h_v22gap);
      ec32gap = bookECorrelatorGap<3,2>("ec32gap",h_v32gap);
      ec42gap = bookECorrelatorGap<4,2>("ec42gap",h_v42gap);
    

      // Integrated, no gap.
      ec22 = bookECorrelator<2,2>("ec22",h_v24);
      ec24 = bookECorrelator<2,4>("ec24",h_v24);
      ec26 = bookECorrelator<2,6>("ec26",h_v26);
      ec28 = bookECorrelator<2,8>("ec28",h_v28);
     
      // pT differential, no gap, three centralities. 
      ec22pT10 = bookECorrelator<2,2>("ec22pT10",h_v24pT10);
      ec24pT10 = bookECorrelator<2,4>("ec24pT10",h_v24pT10);

      ec22pT20 = bookECorrelator<2,2>("ec22pT20",h_v24pT20);
      ec24pT20 = bookECorrelator<2,4>("ec24pT20",h_v24pT20);
    
      ec22pT30 = bookECorrelator<2,2>("ec22pT30",h_v24pT30);
      ec24pT30 = bookECorrelator<2,4>("ec24pT30",h_v24pT30);

      // pT differential, with gap, 30-40% centrality.
      ec22gappT = bookECorrelatorGap<2,2>("ec22gappT",h_v22gappT);
      ec32gappT = bookECorrelatorGap<3,2>("ec32gappT",h_v32gappT);
      ec42gappT = bookECorrelatorGap<4,2>("ec42gappT",h_v42gappT);


      pair<int, int> max = getMaxValues();
      // Declare correlator projections.
      declare(Correlators(cfs, max.first, max.second, h_v22gappT),
        "Correlators");
      declare(Correlators(cfsp, max.first, max.second, h_v22gappT),
        "CorrelatorsPos");
      declare(Correlators(cfsn, max.first, max.second, h_v22gappT),
        "CorrelatorsNeg");

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;

      const double w = event.weight();
      
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
     
      ec22gap->fill(cent, cp, cn, w);
      ec32gap->fill(cent, cp, cn, w); 
      ec42gap->fill(cent, cp, cn, w); 
      
      ec22->fill(cent, c, w);
      ec24->fill(cent, c, w);
      ec26->fill(cent, c, w);
      ec28->fill(cent, c, w);

      if (cent < 10.) ec22pT10->fill(c, w), ec24pT10->fill(c, w);
      else if (cent < 20.) ec22pT20->fill(c, w), ec24pT20->fill(c, w);
      else if (cent < 30.) ec22pT30->fill(c, w), ec24pT30->fill(c, w);
      else if (cent < 40.) {
        ec22gappT->fill(cp, cn, w);
        ec32gappT->fill(cp, cn, w);
        ec42gappT->fill(cp, cn, w);	
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
       // Correlators must be streamed 
       // in order to run reentrant finalize.
       stream();
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
  DECLARE_RIVET_PLUGIN(ALICE_2016_I1419244);


}

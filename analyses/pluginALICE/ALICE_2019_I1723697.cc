// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Correlators.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {


  /// @brief Multiparticle azimuthal correlations pp, pPb, XeXe and PbPb.
  class ALICE_2019_I1723697 : public CumulantAnalysis {
  public:

    /// Constructor
    ALICE_2019_I1723697() : 
      CumulantAnalysis("ALICE_2019_I1723697") {}



    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Find system type.
      string bOpt = getOption<string>("beam","");
      const ParticlePair& beam = beams();
      if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON)
         sysType = pp;
      else if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::LEAD)
        sysType = pPb;
      else if (beam.first.pid() == PID::XENON && beam.second.pid() == PID::XENON)
        sysType = XeXe;
      else if (beam.first.pid() == PID::LEAD && beam.second.pid() == PID::LEAD)
        sysType = PbPb;
      else{
        MSG_WARNING("Suspicious beam. You're probably in reentrant mode. Fetching beam type from option.");
        if (bOpt == "pp") sysType = pp;
        else if (bOpt == "pPb") sysType = pPb;
        else if (bOpt == "XeXe") sysType = XeXe;
        else if (bOpt == "PbPb") sysType = PbPb;
        else MSG_ERROR("Could not decipher beam type. For rivet-merge, set -a ALICE_2019_I1723697:beam=OPT, where opt is pp, pPb, XeXe or PbPb.");        
        MSG_WARNING("Setting beam type to "+bOpt);
      }
      // Sanity check
      if ((sysType == pp && bOpt != "pp") || (sysType == pPb && bOpt != "pPb") ||
        (sysType == XeXe && bOpt != "XeXe") || (sysType == PbPb && bOpt != "PbPb"))
          MSG_WARNING("Beam option and registered beam don't match. Are you sure this is intentional?");

      // Initialise and register projections
      // Declare the trigger projection.
      declare<ALICE::V0AndTrigger>(ALICE::V0AndTrigger(),"V0-AND");
      
      // Centrality projection for high multiplicity trigger in pp.
      if (sysType == pp)
        declareCentrality(ALICE::V0MMultiplicity(), 
	  "ALICE_2015_PPCentrality", "V0M","V0M");

      // The full central charged final state.
      const ChargedFinalState& cfs = ChargedFinalState(Cuts::abseta < 0.8 &&
        Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfs, "CFS");

      // The positive eta side used for rapidity gap = 1.4.
      const ChargedFinalState& cfsp14 = ChargedFinalState(Cuts::eta > 0.5 &&
        Cuts::eta < 0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfsp14, "CFSP14");
      // ..negative ditto.
      const ChargedFinalState& cfsn14 = ChargedFinalState(Cuts::eta < -0.5 &&
        Cuts::eta > -0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfsn14, "CFSN14");
      // The positive eta side used for rapidity gap = 1.0.
      const ChargedFinalState& cfsp10 = ChargedFinalState(Cuts::eta > 0.5 &&
        Cuts::eta < 0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfsp10, "CFSP10");
      // ..negative ditto.
      const ChargedFinalState& cfsn10 = ChargedFinalState(Cuts::eta < -0.5 &&
        Cuts::eta > -0.8 && Cuts::pT > 0.2*GeV && Cuts::pT < 5.0*GeV);
      declare(cfsn10, "CFSN10");

      // Book flow coeff scatters before booking the correlators
      // to have access to bin edges.
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;
      if (sysType == pPb) n1 = 9, n2 = 1;
      else if (sysType == XeXe) n1 = 19, n2 = 1, n3 = 0;
      else if (sysType == PbPb) n1 = 30, n2 = 1, n3 = 1;
      book(h_v22gap, 1 + n1, 1, 1, true);
      book(h_v32gap, 2 + n1, 1, 1, true);
      book(h_v42gap, 3 + n1, 1, 1, true);
      if (sysType != pp)
        book(h_v24, 4 + n1, 1, 1, true);
      book(h_v26, 5 + n1 + n2, 1, 1, true);
      if (sysType == XeXe || sysType == PbPb)
        book(h_v28, 6 + n1 + n2 + n3, 1, 1, true);

      // Book cumulant scatters.
      book(h_c22gap, "c22gap", refData(1 + n1, 1, 1));
      book(h_c32gap, "c32gap", refData(2 + n1, 1, 1));
      book(h_c42gap, "c42gap", refData(3 + n1, 1, 1));
      book(h_c24, "c24", refData(4 + n1, 1, 1));
      book(h_c26, "c26", refData(5 + n1 + n2, 1, 1));
      if (sysType == XeXe || sysType == PbPb)
        book(h_c28, "c28", refData(6 + n1 + n2 + n3, 1, 1));
     
      // Book correlators. First the ungapped ones.
      ec22 = bookECorrelator<2,2>("ec22", refData(4 + n1, 1, 1));
      ec24 = bookECorrelator<2,4>("ec24", refData(4 + n1, 1, 1));
      
      ec622 = bookECorrelator<2,2>("ec622", refData(5 + n1 + n2, 1, 1));
      ec624 = bookECorrelator<2,4>("ec624", refData(5 + n1 + n2, 1, 1));
      ec626 = bookECorrelator<2,6>("ec626", refData(5 + n1 + n2, 1, 1));
      
      // And the ones just valid for XeXe and PbPb.
      if (sysType == XeXe || sysType == PbPb) {
	    ec822 = bookECorrelator<2,2>("ec822", refData(6 + n1 + n2 + n3, 1, 1));
	    ec824 = bookECorrelator<2,4>("ec824", refData(6 + n1 + n2 + n3, 1, 1));
	    ec826 = bookECorrelator<2,6>("ec826", refData(6 + n1 + n2 + n3, 1, 1));
	    ec828 = bookECorrelator<2,8>("ec828", refData(6 + n1 + n2 + n3, 1, 1));
      }
      // ...and the gapped ones.
      ec22gap = bookECorrelatorGap<2,2>("ec22gap", refData(1 + n1, 1, 1));
      ec32gap = bookECorrelatorGap<3,2>("ec32gap", refData(2 + n1, 1, 1));
      ec42gap = bookECorrelatorGap<4,2>("ec42gap", refData(3 + n1, 1, 1));
      // Get the max order of booked correlators for the projections.
      pair<int, int> max = getMaxValues();

      // Declare correlator projections.
      declare(Correlators(cfs, max.first, max.second), "Correlators");
      declare(Correlators(cfsp14, max.first, max.second), "CorrelatorsPos14");
      declare(Correlators(cfsn14, max.first, max.second), "CorrelatorsNeg14");
      declare(Correlators(cfsp10, max.first, max.second), "CorrelatorsPos10");
      declare(Correlators(cfsn10, max.first, max.second), "CorrelatorsNeg10");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;

      // High multiplicity trigger for pp.
      if (sysType == pp && apply<CentralityProjection>(event, "V0M")() > 0.1) 
        vetoEvent;
      // Number of charged particles.
      const double nch = apply<ChargedFinalState>(event, "CFS").size();
      
      // The correlators projections.
      const Correlators& c = applyProjection<Correlators>(event,"Correlators");
      const Correlators& cp10 =
        applyProjection<Correlators>(event,"CorrelatorsPos10");
      const Correlators& cn10 =
        applyProjection<Correlators>(event,"CorrelatorsNeg10");
      const Correlators& cp14 =
        applyProjection<Correlators>(event,"CorrelatorsPos14");
      const Correlators& cn14 =
        applyProjection<Correlators>(event,"CorrelatorsNeg14");

      // Fill correlators.
      ec22->fill(nch, c);
      ec24->fill(nch, c);

      ec622->fill(nch, c);
      ec624->fill(nch, c);
      ec626->fill(nch, c);
      
      if (sysType == XeXe || sysType == PbPb) {
        ec822->fill(nch, c);
        ec824->fill(nch, c);
        ec826->fill(nch, c);
        ec828->fill(nch, c);
      }

      // Fill gapped correlators.
      ec22gap->fill(nch, cp14, cn14);
      ec32gap->fill(nch, cp10, cn10);
      ec42gap->fill(nch, cp10, cn10);

    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Stream correlators to yoda file.
      // Fill cumulant scatters.
      cnTwoInt(h_c22gap, ec22gap);
      cnTwoInt(h_c32gap, ec32gap);
      cnTwoInt(h_c42gap, ec42gap);

      cnFourInt(h_c24, ec22, ec24);
      cnSixInt(h_c26, ec622, ec624, ec626);
      if (sysType == XeXe || sysType == PbPb)
        cnEightInt(h_c28, ec822, ec824, ec826, ec828);

      // Fill flow scatters.
      vnTwoInt(h_v22gap, ec22gap);
      vnTwoInt(h_v32gap, ec32gap);
      vnTwoInt(h_v42gap, ec42gap);
      
      if (sysType != pp)
        vnFourInt(h_v24, ec22, ec24);
      vnSixInt(h_v26, ec622, ec624, ec626);
      if (sysType == XeXe || sysType == PbPb)
        vnEightInt(h_v28, ec822, ec824, ec826, ec828);
    }

    //@}
    // System check enum.
    enum SysType {pp, pPb, XeXe, PbPb};
    SysType sysType;

    /// @name Histograms
    //@{
    // Flow coefficients.
    Scatter2DPtr h_v22gap;
    Scatter2DPtr h_v32gap;
    Scatter2DPtr h_v42gap;
    Scatter2DPtr h_v24;
    Scatter2DPtr h_v26;
    Scatter2DPtr h_v28;

    // Cumulants.
    Scatter2DPtr h_c22gap;
    Scatter2DPtr h_c32gap;
    Scatter2DPtr h_c42gap;
    Scatter2DPtr h_c24;
    Scatter2DPtr h_c26;
    Scatter2DPtr h_c28;

    // Correlators.
    ECorrPtr ec22;
    ECorrPtr ec24;

    ECorrPtr ec622;
    ECorrPtr ec624;
    ECorrPtr ec626;
    
    ECorrPtr ec822;
    ECorrPtr ec824;
    ECorrPtr ec826;
    ECorrPtr ec828;

    // Gapped correlators.
    ECorrPtr ec22gap;
    ECorrPtr ec32gap;
    ECorrPtr ec42gap;

    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2019_I1723697);


}

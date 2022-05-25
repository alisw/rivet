// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class CMS_2018_I1680318 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1680318);


    /// Book histograms and initialise projections before the run
    void init() {

      // Cuts
      MinEnergy     = 5.0; // Particle's energy cut in the forward region [GeV]
      EtaForwardMin = 3.0;
      EtaForwardMax = 5.0;
      EtaCentralCut = 2.4;
      MinParticlePt = 0.5; // [GeV]


      // Initialise and register projections
      const FinalState fsa(Cuts::abseta < EtaForwardMax);
      declare(fsa, "FSA");

      const ChargedFinalState cfs(Cuts::abseta < EtaCentralCut && Cuts::pT > MinParticlePt*GeV);
      declare(cfs, "CFS");

      // Event counters
      book(_num_evts_noCuts,          "num_evts_noCuts");
      book(_num_evts_after_cuts_or,   "num_evts_after_cuts_or");
      book(_num_evts_after_cuts_and,  "num_evts_after_cuts_and");
      book(_num_evts_after_cuts_xor,  "num_evts_after_cuts_xor");
      book(_num_evts_after_cuts_xorm, "num_evts_after_cuts_xorm");
      book(_num_evts_after_cuts_xorp, "num_evts_after_cuts_xorp");

      // Histograms
      book(_hist_dNch_all_dEta_OR,          1,1,1);
      book(_hist_dNch_all_dEta_AND,         1,2,1);
      book(_hist_dNch_all_dEta_XOR,         1,3,1);
      book(_hist_dNch_all_dEta_XORpm,       1,4,1);

      book(_hist_dNch_all_dpt_OR,           2,1,1);
      book(_hist_dNch_all_dpt_AND,          2,2,1);
      book(_hist_dNch_all_dpt_XOR,          2,3,1);

      book(_hist_dNch_leading_dpt_OR,       3,1,1);
      book(_hist_dNch_leading_dpt_AND,      3,2,1);
      book(_hist_dNch_leading_dpt_XOR,      3,3,1);

      book(_hist_integrated_leading_pt_OR,  4,1,1);
      book(_hist_integrated_leading_pt_AND, 4,2,1);
      book(_hist_integrated_leading_pt_XOR, 4,3,1);

      book(_hist_dNev_all_dM_OR,            5,1,1);
      book(_hist_dNev_all_dM_AND,           5,2,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      const FinalState& fsa = apply<FinalState>(event, "FSA");

      bool activity_plus_side  = false,
        activity_minus_side = false;

      for (const Particle& p : fsa.particles()) {
        if ( p.energy() >= MinEnergy ) {
          if ( inRange(p.eta(),      EtaForwardMin,      EtaForwardMax) ) activity_plus_side  = true;
          if ( inRange(p.eta(), -1.0*EtaForwardMax, -1.0*EtaForwardMin) ) activity_minus_side = true;
        }

        // If activity already found in both sides,
        // then there is no point in keep going the loop
        if (activity_plus_side && activity_minus_side) break;
      }

      // Event selections
      const bool cutsor   = ( activity_plus_side ||  activity_minus_side);
      const bool cutsand  = ( activity_plus_side &&  activity_minus_side);
      const bool cutsxor  = ((activity_plus_side && !activity_minus_side) || (!activity_plus_side && activity_minus_side));
      const bool cutsxorm = (!activity_plus_side &&  activity_minus_side);
      const bool cutsxorp = ( activity_plus_side && !activity_minus_side);

      _num_evts_noCuts->fill();
      if ( charged.size() >= 1 ) {
        if (cutsor)   _num_evts_after_cuts_or   ->fill();
        if (cutsand)  _num_evts_after_cuts_and  ->fill();
        if (cutsxor)  _num_evts_after_cuts_xor  ->fill();
        if (cutsxorm) _num_evts_after_cuts_xorm ->fill();
        if (cutsxorp) _num_evts_after_cuts_xorp ->fill();
      }

      // Loop over charged particles
      double leading_pt = 0;
      for (const Particle& p : charged.particles()) {
        // Find the leading-pt particle of the event
        if (p.pT() > leading_pt) leading_pt = p.pT();

        // Filling histograms
        if (cutsor)   _hist_dNch_all_dEta_OR    -> fill(p.eta());
        if (cutsand)  _hist_dNch_all_dEta_AND   -> fill(p.eta());
        if (cutsxor)  _hist_dNch_all_dEta_XOR   -> fill(p.eta());

        //Average xorm & xorp
        if (cutsxorm) _hist_dNch_all_dEta_XORpm -> fill(p.eta());
        if (cutsxorp) _hist_dNch_all_dEta_XORpm -> fill(-1.0*p.eta());

        if (cutsor)   _hist_dNch_all_dpt_OR     -> fill(p.pT());
        if (cutsand)  _hist_dNch_all_dpt_AND    -> fill(p.pT());
        if (cutsxor)  _hist_dNch_all_dpt_XOR    -> fill(p.pT());
      }

      // Filling multiplicity histograms
      if ( charged.size() >= 1 ) {
        if (cutsor)  _hist_dNev_all_dM_OR  -> fill(charged.size());
        if (cutsand) _hist_dNev_all_dM_AND -> fill(charged.size());
      }

      // Filling leading-pt histograms
      if (cutsor)  _hist_dNch_leading_dpt_OR  -> fill(leading_pt);
      if (cutsand) _hist_dNch_leading_dpt_AND -> fill(leading_pt);
      if (cutsxor) _hist_dNch_leading_dpt_XOR -> fill(leading_pt);

      // Integrating leading-pt histograms
      for (size_t i = 0 ; i < _hist_integrated_leading_pt_OR->numBins() ; ++i) {
        double binlimitlow_t = _hist_integrated_leading_pt_OR->bin(i).xMin();
        double weightbw_t    = _hist_integrated_leading_pt_OR->bin(i).xWidth();
        double xbin_t        = _hist_integrated_leading_pt_OR->bin(i).xMid();
        if (leading_pt > binlimitlow_t) {
          if (cutsor)  _hist_integrated_leading_pt_OR  -> fill(xbin_t, weightbw_t);
          if (cutsand) _hist_integrated_leading_pt_AND -> fill(xbin_t, weightbw_t);
          if (cutsxor) _hist_integrated_leading_pt_XOR -> fill(xbin_t, weightbw_t);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_INFO("Number of selected events: "                  << endl
               << "\t All       = " << _num_evts_noCuts->val()          << endl
               << "\t Inelastic = " << _num_evts_after_cuts_or->val()   << endl
               << "\t NSD       = " << _num_evts_after_cuts_and->val()  << endl
               << "\t Xor       = " << _num_evts_after_cuts_xor->val()  << endl
               << "\t Xorm      = " << _num_evts_after_cuts_xorm->val() << endl
               << "\t Xorp      = " << _num_evts_after_cuts_xorp->val());

      scale(_hist_dNch_all_dEta_OR,    1./ *_num_evts_after_cuts_or);
      scale(_hist_dNch_all_dEta_AND,   1./ *_num_evts_after_cuts_and);
      scale(_hist_dNch_all_dEta_XOR,   1./ *_num_evts_after_cuts_xor);
      scale(_hist_dNch_all_dEta_XORpm, 1./ (*_num_evts_after_cuts_xorm + *_num_evts_after_cuts_xorp));

      scale(_hist_dNch_all_dpt_OR,   1./ *_num_evts_after_cuts_or);
      scale(_hist_dNch_all_dpt_AND,  1./ *_num_evts_after_cuts_and);
      scale(_hist_dNch_all_dpt_XOR,  1./ *_num_evts_after_cuts_xor);

      scale(_hist_dNch_leading_dpt_OR,   1./ *_num_evts_after_cuts_or);
      scale(_hist_dNch_leading_dpt_AND,  1./ *_num_evts_after_cuts_and);
      scale(_hist_dNch_leading_dpt_XOR,  1./ *_num_evts_after_cuts_xor);

      scale(_hist_integrated_leading_pt_OR,   1./ *_num_evts_after_cuts_or);
      scale(_hist_integrated_leading_pt_AND,  1./ *_num_evts_after_cuts_and);
      scale(_hist_integrated_leading_pt_XOR,  1./ *_num_evts_after_cuts_xor);

      scale(_hist_dNev_all_dM_OR,   1./ *_num_evts_after_cuts_or);
      scale(_hist_dNev_all_dM_AND,  1./ *_num_evts_after_cuts_and);
    }


  private:

    // Cuts
    double MinEnergy, EtaForwardMin, EtaForwardMax, EtaCentralCut, MinParticlePt;

    // Counters
    CounterPtr _num_evts_noCuts,
                  _num_evts_after_cuts_and,
                  _num_evts_after_cuts_or,
                  _num_evts_after_cuts_xor,
                  _num_evts_after_cuts_xorp,
                  _num_evts_after_cuts_xorm;

    // Histograms
    Histo1DPtr
    _hist_dNch_all_dEta_AND,
                  _hist_dNch_all_dEta_OR,
                  _hist_dNch_all_dEta_XOR,
                  _hist_dNch_all_dEta_XORpm;
    Histo1DPtr
    _hist_dNch_all_dpt_AND,
                  _hist_dNch_all_dpt_OR,
                  _hist_dNch_all_dpt_XOR;
    Histo1DPtr
    _hist_dNch_leading_dpt_AND,
                  _hist_dNch_leading_dpt_OR,
                  _hist_dNch_leading_dpt_XOR;
    Histo1DPtr
    _hist_integrated_leading_pt_AND,
                  _hist_integrated_leading_pt_OR,
                  _hist_integrated_leading_pt_XOR;
    Histo1DPtr
    _hist_dNev_all_dM_AND,
                  _hist_dNev_all_dM_OR;

  };


  RIVET_DECLARE_PLUGIN(CMS_2018_I1680318);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  /// @brief BRAHMS Centrality projection.
  class BRAHMSCentrality : public SingleValueProjection {
  public:
    // Constructor
    BRAHMSCentrality() : SingleValueProjection() {
      // Using here the BRAHMS reaction centrality from eg. 1602.01183, which
      // might not be correct.
      declare(ChargedFinalState(Cuts::pT > 0.1*GeV && Cuts::abseta < 2.2),
        "ChargedFinalState");
    }
    // Destructor
    virtual ~BRAHMSCentrality() {}

    // Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(BRAHMSCentrality);

  protected:
    // Do the projection. Count the number of charged particles in
    // the specified range.
    virtual void project(const Event& e) {
      clear();
      set(apply<ChargedFinalState>
        (e, "ChargedFinalState").particles().size());
    }

    // Compare to another projection.
    virtual int compare(const Projection& p) const {
      // This projection is only used for the analysis below.
      return UNDEFINED;
    }

  };

  /// @brief Brahms centrality calibration analysis based on the
  //  BrahmsCentrality projection. No data is given for this
  //  analysis, so one MUST do a calibration run.
  class BRAHMS_2004_CENTRALITY : public Analysis {
  public:
    // Constructor
    BRAHMS_2004_CENTRALITY() : Analysis("BRAHMS_2004_CENTRALITY") {}

    // Initialize the analysis
    void init() {
       declare(BRAHMSCentrality(),"Centrality");
       declare(ImpactParameterProjection(), "IMP");
       
       // The central multiplicity.
       mult = bookHisto1D("mult",450,0,4500);
       
       // Safeguard against filling preloaded histograms.
       done = (mult->numEntries() > 0);

       // The impact parameter.
       imp = bookHisto1D("mult_IMP",100,0,20);
    }

    // Analyse a single event
    void analyze(const Event& event) {
      if (done) return;
      // Fill impact parameter.
      imp->fill(apply<SingleValueProjection>(event,"IMP")(), event.weight());
      // Fill multiplicity.
      mult->fill(apply<SingleValueProjection>(event,"Centrality")(), event.weight());
    }

    // Finalize the analysis
    void finalize() {
      // Normalize the distributions, safeguarding against 
      // yoda normalization error.
      if(mult->numEntries() > 0) mult->normalize();
      if(imp->numEntries() > 0) imp->normalize();
    
    }

  private:
    // Histograms.
    Histo1DPtr mult;
    Histo1DPtr imp;
    // Flag to test if we have preloaded histograms.
    bool done;
  
  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BRAHMS_2004_CENTRALITY);

  /// @brief Brahms pT spectra for id particles (pi+, pi-, K+, K-)
  //  in small bins of rapidity, 5% central collisions. 
  //  System: AuAu @ 200GeV/nn.
  class BRAHMS_2004_I647076 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BRAHMS_2004_I647076);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Centrality Projection.
      declareCentrality(BRAHMSCentrality(), "BRAHMS_2004_CENTRALITY","mult","BCEN");
      // TODO: Feed down correction is unclear.
      declare(FinalState(Cuts::rap < 4 && Cuts::rap > -0.1 && Cuts::pT > 100*MeV), "FS");
      // The measured rapidity intervals for pions.
      rapIntervalsPi = {{-0.1,0.},{0.,0.1},{0.4,0.6},{0.6,0.8},{0.8,1.0},
        {1.0,1.2},{1.2,1.4},{2.1,2.3},{2.4,2.6},{3.0,3.1},{3.1,3.2},{3.2,3.3},
        {3.3,3.4},{3.4,3.66}};
      // The measured rapidity intervals for kaons.
      rapIntervalsK = {{-0.1,0.},{0.,0.1},{0.4,0.6},{0.6,0.8},{0.8,1.0},
        {1.0,1.2},{2.0,2.2},{2.3,2.5},{2.9,3.0},{3.0,3.1},{3.1,3.2},{3.2,3.4}};
      // Book histograms
      for (int i = 1, N = rapIntervalsPi.size(); i <= N; ++i) {
        piPlus.push_back(bookHisto1D(1, 1, i));
        piMinus.push_back(bookHisto1D(1, 1, 14 + i));
      }
      for (int i = 1, N = rapIntervalsK.size(); i <= N; ++i) {
        kPlus.push_back(bookHisto1D(2, 1, i));
        kMinus.push_back(bookHisto1D(2, 1, 12 + i));
      }
      // Counter for accepted sum of weights (centrality cut).
      centSow = bookCounter("centSow");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double w = event.weight();
      // Reject all non-central events. The paper does not speak of 
      // any other event trigger, which in any case should matter
      // little for central events.
      if(apply<CentralityProjection>(event,"BCEN")() > 5.0) return;
      // Keep track of sum of weights.
      centSow->fill(w);
      const FinalState& fs = apply<FinalState>(event,"FS");
      // Loop over particles.
      for (const auto& p : fs.particles()) {
        const double y = p.rapidity();
	const double pT = p.pT();
	const int id = p.pid();
	// First pions.
	if (abs(id) == 211) {
          // Protect against decaying K0S and Lambda
	  if (p.hasAncestor(310) || p.hasAncestor(-310) || 
	    p.hasAncestor(3122) || p.hasAncestor(3122)) continue;
	  for (int i = 0, N = rapIntervalsPi.size(); i < N; ++i) {
	    if (y > rapIntervalsPi[i].first && y <= rapIntervalsPi[i].second) {
	      const double dy = rapIntervalsPi[i].second - rapIntervalsPi[i].first;
	      const double nWeight = w / ( 2.*M_PI*pT*dy);
	      if (id == 211) piPlus[i]->fill(pT, nWeight);
	      else piMinus[i]->fill(pT, nWeight);
	      break;
	    }
	  } 
	}
	// Then kaons.
	else if (abs(id) == 321) {
	  for (int i = 0, N = rapIntervalsK.size(); i < N; ++i) {
	    if (y > rapIntervalsK[i].first && y <= rapIntervalsK[i].second) {
	      const double dy = rapIntervalsK[i].second - rapIntervalsK[i].first;
	      const double nWeight = w / ( 2.*M_PI*pT*dy);
	      if (id == 321) kPlus[i]->fill(pT, nWeight);
	      else kMinus[i]->fill(pT, nWeight);
	      break;
	    }
	  } 
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Normalize all histograms to per-event yields.
      for (int i = 0, N = rapIntervalsPi.size(); i < N; ++i) {
        piPlus[i]->scaleW(1./centSow->sumW());
        piMinus[i]->scaleW(1./centSow->sumW());
      }
      for (int i = 0, N = rapIntervalsK.size(); i < N; ++i) {
        kPlus[i]->scaleW(1./centSow->sumW());
        kMinus[i]->scaleW(1./centSow->sumW());
      }

    }

    //@}

    // The rapidity intervals.
    vector<pair<double, double> > rapIntervalsPi;
    vector<pair<double, double> > rapIntervalsK;

    /// @name Histograms
    //@{
    vector<Histo1DPtr> piPlus;
    vector<Histo1DPtr> piMinus;
    vector<Histo1DPtr> kPlus;
    vector<Histo1DPtr> kMinus;
    CounterPtr centSow;
    //@}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BRAHMS_2004_I647076);
}

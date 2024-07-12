// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#ifdef RIVET_ENABLE_HEPMC_3    
#include "HepMC3/WriterAscii.h"
#endif

namespace Rivet {


  /// @brief Just measures a few observables as a demo and 
  /// writes to a file all the events with thrust < 0.1
  class EXAMPLE_HEPMC_FILTER : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_HEPMC_FILTER);

    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Projections
      const FinalState fs(Cuts::abseta < 2.5);
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");

      #ifdef RIVET_ENABLE_HEPMC_3   
      _writer = std::make_shared<HepMC3::WriterAscii>("EXAMPLE_HEPMC_FILTER.hepmc3");
      #else
      _writer = std::make_shared<HepMC::IO_GenEvent>("EXAMPLE_HEPMC_FILTER.hepmc2", std::ios::out);
      #endif
    }


    /// Do the analysis
    void analyze(const Event& event) {

      const Jets& jets = apply<FastJets>(event, "Jets").jets(Cuts::pT > 20*GeV);
      const size_t num_b_jets = count(jets, hasBTag(Cuts::pT > 500*MeV));

      if (num_b_jets > 0) {
        #ifdef RIVET_ENABLE_HEPMC_3  
       if (!_writer->run_info()) _writer-> set_run_info(event.genEvent()->run_info());
       _writer->write_event(*event.genEvent());		  
        #else
       _writer->write_event(event.genEvent());		  
        #endif
      }		  
    }


    //@}

    /// @name Output handler
    //@{

    #ifdef RIVET_ENABLE_HEPMC_3    
    std::shared_ptr<RivetHepMC::WriterAscii> _writer;
    #else
    std::shared_ptr<RivetHepMC::IO_GenEvent> _writer;
    #endif

    //@}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(EXAMPLE_HEPMC_FILTER);

}

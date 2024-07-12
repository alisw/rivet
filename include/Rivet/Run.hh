// -*- C++ -*-
#ifndef RIVET_Run_HH
#define RIVET_Run_HH

#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {


  // Forward declaration
  class AnalysisHandler;


  /// @brief Interface to handle a run of events read from a HepMC stream or file.
  class Run {
  public:

    /// Standard constructor.
    Run(AnalysisHandler& ah);

    /// Destructor
    ~Run();


    /// @name Set run properties
    /// @{

    /// Get the cross-section for this run.
    Run& setCrossSection(double xs);

    /// Declare whether to list available analyses
    Run& setListAnalyses(bool dolist);

    /// @}


    /// @name File processing stages
    /// @{

    /// Set up HepMC file readers (using the appropriate file weight for the first file)
    bool init(const std::string& evtfile, double weight=1.0);

    /// Open a HepMC GenEvent file (using the appropriate file weight for the first file)
    bool openFile(const std::string& evtfile, double weight=1.0);

    /// Read the next HepMC event
    bool readEvent();

    /// Read the next HepMC event only to skip it
    //bool skipEvent();

    /// Return the number of (collapsed) events read in from HepMC,
    /// including current partial event in case of sub-events
    size_t numEvents() const { return _evtcount; }

    /// Handle next event
    bool processEvent();

    /// Close up HepMC I/O
    bool finalize();

    /// @}


  private:

    /// Get a Log object
    Log& getLog() const;

    /// AnalysisHandler object
    AnalysisHandler& _ah;

    /// @name Run variables obtained from events or command line
    /// @{

    /// @brief An extra event weight scaling per event file.
    /// Useful for e.g. AlpGen n-parton event file combination.
    double _fileweight = 1.0;

    /// Cross-section from command line.
    double _xs = NAN;

    /// Number of (collapsed) events read from file so far
    size_t _evtcount = 0;

    /// Current event number to keep track of sub-events
    int _evtnumber = -1;

    /// @}


    /// Flag to show list of analyses
    bool _listAnalyses = false;


    /// @name HepMC I/O members
    /// @{

    /// Current event
    std::shared_ptr<GenEvent> _evt;

    /// Output stream for HepMC writer
    std::shared_ptr<std::istream> _istr;

    /// HepMC reader
    std::shared_ptr<HepMC_IO_type> _hepmcReader;

    /// @}

  };


}

#endif

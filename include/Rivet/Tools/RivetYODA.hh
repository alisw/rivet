#ifndef RIVET_RIVETYODA_HH
#define RIVET_RIVETYODA_HH

/// @author Andy Buckley
/// @date   2009-01-30
/// @author David Grellscheid
/// @date   2011-07-18

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetBoost.hh"

#include "YODA/AnalysisObject.h"
#include "YODA/WriterYODA.h"
#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Point2D.h"
#include <map>

namespace Rivet {


  typedef shared_ptr<YODA::AnalysisObject> AnalysisObjectPtr;
  typedef shared_ptr<YODA::Counter> CounterPtr;
  typedef shared_ptr<YODA::Histo1D> Histo1DPtr;
  typedef shared_ptr<YODA::Histo2D> Histo2DPtr;
  typedef shared_ptr<YODA::Profile1D> Profile1DPtr;
  typedef shared_ptr<YODA::Profile2D> Profile2DPtr;
  typedef shared_ptr<YODA::Scatter2D> Scatter2DPtr;
  typedef shared_ptr<YODA::Scatter3D> Scatter3DPtr;

  using YODA::WriterYODA;
  using YODA::Counter;
  using YODA::Histo1D;
  using YODA::HistoBin1D;
  using YODA::Histo2D;
  using YODA::HistoBin2D;
  using YODA::Profile1D;
  using YODA::ProfileBin1D;
  using YODA::Profile2D;
  using YODA::ProfileBin2D;
  using YODA::Scatter1D;
  using YODA::Point1D;
  using YODA::Scatter2D;
  using YODA::Point2D;
  using YODA::Scatter3D;
  using YODA::Point3D;

  /// Function to get a map of all the refdata in a paper with the
  /// given @a papername.
  map<string, Scatter2DPtr> getRefData(const string& papername);

  /// @todo Also provide a Scatter3D getRefData() version?

  /// Get the file system path to the reference file for this paper.
  string getDatafilePath(const string& papername);

  /// Return the integral over the histogram bins
  /// @deprecated Prefer to directly use the histo's integral() method.
  DEPRECATED("Prefer to directly use the histo's integral() method.")
  inline double integral(Histo1DPtr histo) {
    return histo->integral();
  }


}

#endif

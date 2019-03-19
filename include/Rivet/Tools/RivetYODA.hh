#ifndef RIVET_RIVETYODA_HH
#define RIVET_RIVETYODA_HH

#include "Rivet/Config/RivetCommon.hh"
#include "YODA/AnalysisObject.h"
#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

namespace Rivet {

  typedef std::shared_ptr<YODA::AnalysisObject> AnalysisObjectPtr;
  typedef std::shared_ptr<YODA::Counter> CounterPtr;
  typedef std::shared_ptr<YODA::Histo1D> Histo1DPtr;
  typedef std::shared_ptr<YODA::Histo2D> Histo2DPtr;
  typedef std::shared_ptr<YODA::Profile1D> Profile1DPtr;
  typedef std::shared_ptr<YODA::Profile2D> Profile2DPtr;
  typedef std::shared_ptr<YODA::Scatter1D> Scatter1DPtr;
  typedef std::shared_ptr<YODA::Scatter2D> Scatter2DPtr;
  typedef std::shared_ptr<YODA::Scatter3D> Scatter3DPtr;

  using YODA::AnalysisObject;
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


  /// Function to get a map of all the refdata in a paper with the given @a papername.
  map<string, AnalysisObjectPtr> getRefData(const string& papername);

  /// Get the file system path to the reference file for this paper.
  string getDatafilePath(const string& papername);


  /// Traits class to access the type of the AnalysisObject in the
  /// reference files.
  template<typename T> struct ReferenceTraits {};
  template<> struct ReferenceTraits<Counter> { typedef Counter RefT; };
  template<> struct ReferenceTraits<Scatter1D> { typedef Scatter1D RefT; };
  template<> struct ReferenceTraits<Histo1D> { typedef Scatter2D RefT; };
  template<> struct ReferenceTraits<Profile1D> { typedef Scatter2D RefT; };
  template<> struct ReferenceTraits<Scatter2D> { typedef Scatter2D RefT; };
  template<> struct ReferenceTraits<Histo2D> { typedef Scatter3D RefT; };
  template<> struct ReferenceTraits<Profile2D> { typedef Scatter3D RefT; };
  template<> struct ReferenceTraits<Scatter3D> { typedef Scatter3D RefT; };


  /// If @a dst and @a src both are of same subclass T, copy the
  /// contents of @a src into @a dst and return true. Otherwise return
  /// false.
  template <typename T>
  inline bool aocopy(AnalysisObjectPtr src, AnalysisObjectPtr dst) {
    shared_ptr<T> tsrc = dynamic_pointer_cast<T>(src);
    if ( !tsrc ) return false;
    shared_ptr<T> tdst = dynamic_pointer_cast<T>(dst);
    if ( !tdst ) return false;
    *tdst = *tsrc;
    return true;
  }

  /// If @a dst and @a src both are of same subclass T, add the
  /// contents of @a src into @a dst and return true. Otherwise return
  /// false.
  template <typename T>
  inline bool aoadd(AnalysisObjectPtr dst, AnalysisObjectPtr src, double scale) {
    shared_ptr<T> tsrc = dynamic_pointer_cast<T>(src);
    if ( !tsrc ) return false;
    shared_ptr<T> tdst = dynamic_pointer_cast<T>(dst);
    if ( !tdst ) return false;
    tsrc->scaleW(scale);
    *tdst += *tsrc;
    return true;
  }

  /// If @a dst is the same subclass as @a src, copy the contents of @a
  /// src into @a dst and return true. Otherwise return false.
  bool copyao(AnalysisObjectPtr src, AnalysisObjectPtr dst);

  /// If @a dst is the same subclass as @a src, scale the contents of
  /// @a src with @a scale and add it to @a dst and return true. Otherwise
  /// return false.
  bool addaos(AnalysisObjectPtr dst, AnalysisObjectPtr src, double scale);

  /// Check if two analysis objects have the same binning or, if not
  /// binned, are in other ways compatible.
  // inline bool bookingCompatible(CounterPtr, CounterPtr) {
  //   return true;
  // }
  template <typename TPtr>
  inline bool bookingCompatible(TPtr a, TPtr b) {
    return a->sameBinning(*b);
  }
  inline bool bookingCompatible(CounterPtr, CounterPtr) {
    return true;
  }
  inline bool bookingCompatible(Scatter1DPtr a, Scatter1DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(Scatter2DPtr a, Scatter2DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(Scatter3DPtr a, Scatter3DPtr b) {
    return a->numPoints() == b->numPoints();
  }

}

#endif

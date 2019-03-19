#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "YODA/ReaderYODA.h"
#include "YODA/ReaderAIDA.h"

using namespace std;

namespace Rivet {


  string getDatafilePath(const string& papername) {
    /// Try to find YODA otherwise fall back to try AIDA
    const string path1 = findAnalysisRefFile(papername + ".yoda");
    if (!path1.empty()) return path1;
    const string path2 = findAnalysisRefFile(papername + ".aida");
    if (!path2.empty()) return path2;
    throw Rivet::Error("Couldn't find ref data file '" + papername + ".yoda" +
                       " in data path, '" + getRivetDataPath() + "', or '.'");
  }


  map<string, AnalysisObjectPtr> getRefData(const string& papername) {
    const string datafile = getDatafilePath(papername);

    // Make an appropriate data file reader and read the data objects
    /// @todo Remove AIDA support some day...
    YODA::Reader& reader = (datafile.find(".yoda") != string::npos) ?   \
      YODA::ReaderYODA::create() : YODA::ReaderAIDA::create();
    vector<YODA::AnalysisObject *> aovec;
    reader.read(datafile, aovec);

    // Return value, to be populated
    map<string, AnalysisObjectPtr> rtn;
    foreach ( YODA::AnalysisObject* ao, aovec ) {
      AnalysisObjectPtr refdata(ao);
      if (!refdata) continue;
      const string plotpath = refdata->path();
      // Split path at "/" and only return the last field, i.e. the histogram ID
      const size_t slashpos = plotpath.rfind("/");
      const string plotname = (slashpos+1 < plotpath.size()) ? plotpath.substr(slashpos+1) : "";
      rtn[plotname] = refdata;
    }
    return rtn;
  }

  bool copyao(AnalysisObjectPtr src, AnalysisObjectPtr dst) {
    for (const std::string& a : src->annotations())
      dst->setAnnotation(a, src->annotation(a));
    if ( aocopy<Counter>(src,dst) ) return true;
    if ( aocopy<Histo1D>(src,dst) ) return true;
    if ( aocopy<Histo2D>(src,dst) ) return true;
    if ( aocopy<Profile1D>(src,dst) ) return true;
    if ( aocopy<Profile2D>(src,dst) ) return true;
    if ( aocopy<Scatter1D>(src,dst) ) return true;
    if ( aocopy<Scatter2D>(src,dst) ) return true;
    if ( aocopy<Scatter3D>(src,dst) ) return true;
    return false;
  }

  bool addaos(AnalysisObjectPtr dst, AnalysisObjectPtr src, double scale) {
    if ( aoadd<Counter>(dst,src,scale) ) return true;
    if ( aoadd<Histo1D>(dst,src,scale) ) return true;
    if ( aoadd<Histo2D>(dst,src,scale) ) return true;
    if ( aoadd<Profile1D>(dst,src,scale) ) return true;
    if ( aoadd<Profile2D>(dst,src,scale) ) return true;
    return false;
  }


}

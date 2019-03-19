#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/Percentile.hh"
#include "Rivet/Analysis.hh"

using namespace std;

namespace Rivet {

void PercentileBase::selectBins(const Event & ev) {
  const CentralityProjection & proj =
    _ana->apply<CentralityProjection>(ev, _projName);
  _activeBins.clear();
  const int nCent = _cent.size();
  for (int iCent = 0; iCent < nCent; ++iCent) {
    if ( inRange(proj(), _cent[iCent]) )
      _activeBins.push_back(iCent);
  }
}

}

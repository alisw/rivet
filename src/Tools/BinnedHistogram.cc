// -*- C++ -*-
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Analysis.hh"

namespace Rivet {


  const BinnedHistogram& BinnedHistogram::add(double binMin, double binMax, Histo1DPtr histo) {
    if (binMin > binMax) throw RangeError("Cannot add a binned histogram where the lower bin edge is above the upper edge");
    _histosByUpperBound[binMax] = histo;
    _histosByLowerBound[binMin] = histo;
    bool found = false;
    for (Histo1DPtr hist : _histos) {
      if (hist == histo) {
        found = true;
        break;
      }
    }

    if (!found){
      _histos.push_back(histo);
      _binWidths[histo] = binMax-binMin;
    }

    return *this;
  }


  const Histo1DPtr BinnedHistogram::histo(double binval) const {
    // Check that the bin is not out of range
    auto histIt1 = _histosByUpperBound.upper_bound(binval);
    if (histIt1 == _histosByUpperBound.end()) throw RangeError("BinnedHistogram: no bin found");

    Histo1DPtr histo = histIt1->second;

    // No need to check going beyond the upper bound if we already passed above
    // (given that upper bound > lower bound is checked)
    // Check it is not before the start of the map
    auto histIt2 = _histosByLowerBound.lower_bound(binval);
    if (histIt2 == _histosByLowerBound.begin()) throw RangeError("BinnedHistogram: no bin found");

    // By-lower-bound actually gives us the iterator one above the nearest element,
    // so decrement it. This is safe because we already checked we're not at the start!
    --histIt2;
    if (histo != histIt2->second) throw RangeError("BinnedHistogram: no bin found");

    return histo;
  }


  Histo1DPtr BinnedHistogram::histo(double binval) {
    return static_cast<const BinnedHistogram&>(*this).histo(binval); //< trick to avoid duplication
  }


  void BinnedHistogram::fill(double binval, double val, double weight) {
    try {
      Histo1DPtr h = histo(binval);
      h->fill(val, weight);
      //return h;
    } catch (RangeError& e) {} //< no bin found: do nothing
  }


  void BinnedHistogram::scale(double scale, Analysis* ana) {
    for (Histo1DPtr hist : histos()) {
      ana->scale(hist, scale/_binWidths[hist]);
    }
  }


}

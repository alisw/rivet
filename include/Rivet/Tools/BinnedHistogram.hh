// -*- C++ -*-
#ifndef RIVET_BINNEDHISTOGRAM_HH
#define RIVET_BINNEDHISTOGRAM_HH

/// @todo BinnedHistogram needs to have a list of interbnal members first which then get booked
/// by the analysis. Booking a temporary, and then adding into BinnedHisto is not possible


#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {

  class Analysis;


  /// @brief A set of booked Histo1DPtr, each in a bin of a second variable.
  ///
  /// BinnedHistogram contains a series of histograms of the same quantity
  /// each in a different region of a second quantity.  For example, a
  /// BinnedHistogram may contain histograms of the cross-section differential
  /// in \f$ p_T \f$ in different \f$ \eta \f$  regions.
  class BinnedHistogram {
  public:

    /// Create a new empty BinnedHistogram
    BinnedHistogram() = default;

    /// Create a new BinnedHistogram with the given bin edges and contents
    BinnedHistogram(const vector<double>& edges, const vector<Histo1DPtr>& histos) {
      assert(edges.size() == histos.size()+1);
      for (size_t i = 0; i < histos.size(); ++i)
        add(edges[i], edges[i+1], histos[i]);
    }

    /// @todo Can we have an "emplace constructor", passing tuples of bookHisto1D args?


    ///  Add a histogram in the @c T bin between @a binMin and @a binMax
    const BinnedHistogram & add(double binMin, double binMax, Histo1DPtr histo);


    /// Fill the histogram in the same bin as @a binval with value @a val and weight @a weight
    void fill(double binval, double val, double weight = 1.0);


    /// @brief Get the histogram in the same bin as @a binval (const)
    /// @note Throws a RangeError if @a binval doesn't fall in a declared bin
    const Histo1DPtr histo(double binval) const;
    /// @brief Get the histogram in the same bin as @a binval
    /// @note Throws a RangeError if @a binval doesn't fall in a declared bin
    Histo1DPtr histo(double binval);

    /// Get the contained histograms (const)
    const vector<Histo1DPtr>& histos() const { return _histos; }
    /// Get the contained histograms
    vector<Histo1DPtr>& histos() { return _histos; }



    /// Scale histograms taking into account its "external" binwidth, i.e. by scale/binWidth
    /// @note The Analysis pointer is passed in order to call the analysis' scale(h) method: can we avoid that?
    void scale(double scale, Analysis* ana);


  private:

    map<double, Histo1DPtr> _histosByUpperBound, _histosByLowerBound;
    vector<Histo1DPtr> _histos;
    map<Histo1DPtr, double> _binWidths;

  };


}

#endif

// -*- C++ -*-
#ifndef RIVET_Analysis_HH
#define RIVET_Analysis_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleUtils.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Tools/RivetMT2.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/Percentile.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include <tuple>


/// @def vetoEvent
/// Preprocessor define for vetoing events, including the log message and return.
#define vetoEvent                                                       \
  do { MSG_DEBUG("Vetoing event on line " << __LINE__ << " of " << __FILE__); return; } while(0)


namespace Rivet {


  // Convenience for analysis writers
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::tuple;
  using std::stringstream;
  using std::swap;
  using std::numeric_limits;


  // Forward declaration
  class AnalysisHandler;


  /// @brief This is the base class of all analysis classes in Rivet.
  ///
  /// There are
  /// three virtual functions which should be implemented in base classes:
  ///
  /// void init() is called by Rivet before a run is started. Here the
  /// analysis class should book necessary histograms. The needed
  /// projections should probably rather be constructed in the
  /// constructor.
  ///
  /// void analyze(const Event&) is called once for each event. Here the
  /// analysis class should apply the necessary Projections and fill the
  /// histograms.
  ///
  /// void finalize() is called after a run is finished. Here the analysis
  /// class should do whatever manipulations are necessary on the
  /// histograms. Writing the histograms to a file is, however, done by
  /// the Rivet class.
  class Analysis : public ProjectionApplier {
  public:

    /// The AnalysisHandler is a friend.
    friend class AnalysisHandler;


    /// Constructor
    Analysis(const std::string& name);

    /// The destructor
    virtual ~Analysis() {}

    /// The assignment operator is private and mustdeleted, so it can never be called.
    Analysis& operator=(const Analysis&) = delete;


  public:

    /// @defgroup analysis_main Main analysis methods
    /// @{

    /// Initialize this analysis object. A concrete class should here
    /// book all necessary histograms. An overridden function must make
    /// sure it first calls the base class function.
    virtual void init() { }

    /// Analyze one event. A concrete class should here apply the
    /// necessary projections on the \a event and fill the relevant
    /// histograms. An overridden function must make sure it first calls
    /// the base class function.
    virtual void analyze(const Event& event) = 0;

    /// Finalize this analysis object. A concrete class should here make
    /// all necessary operations on the histograms. Writing the
    /// histograms to a file is, however, done by the Rivet class. An
    /// overridden function must make sure it first calls the base class
    /// function.
    virtual void finalize() { }

    /// @}


  public:

    /// @defgroup analysis_meta Metadata
    /// Metadata is used for querying from the command line and also for
    /// building web pages and the analysis pages in the Rivet manual.
    /// @{

    /// Get the actual AnalysisInfo object in which all this metadata is stored.
    const AnalysisInfo& info() const {
      assert(_info && "No AnalysisInfo object :O");
      return *_info;
    }

    /// @brief Get the name of the analysis.
    ///
    /// By default this is computed by combining the results of the
    /// experiment, year and Spires ID metadata methods and you should
    /// only override it if there's a good reason why those won't
    /// work. If options has been set for this instance, a
    /// corresponding string is appended at the end.
    virtual std::string name() const {
      return  ( (info().name().empty()) ? _defaultname : info().name() ) + _optstring;
    }

    /// Get name of reference data file, which could be different from plugin name
    virtual std::string getRefDataName() const {
      return (info().getRefDataName().empty()) ? _defaultname : info().getRefDataName();
    }

    /// Set name of reference data file, which could be different from plugin name
    virtual void setRefDataName(const std::string& ref_data="") {
      info().setRefDataName(!ref_data.empty() ? ref_data : name());
    }

    /// Get the Inspire ID code for this analysis.
    virtual std::string inspireId() const {
      return info().inspireId();
    }

    /// Get the SPIRES ID code for this analysis (~deprecated).
    virtual std::string spiresId() const {
      return info().spiresId();
    }

    /// @brief Names & emails of paper/analysis authors.
    ///
    /// Names and email of authors in 'NAME \<EMAIL\>' format. The first
    /// name in the list should be the primary contact person.
    virtual std::vector<std::string> authors() const {
      return info().authors();
    }

    /// @brief Get a short description of the analysis.
    ///
    /// Short (one sentence) description used as an index entry.
    /// Use @a description() to provide full descriptive paragraphs
    /// of analysis details.
    virtual std::string summary() const {
      return info().summary();
    }

    /// @brief Get a full description of the analysis.
    ///
    /// Full textual description of this analysis, what it is useful for,
    /// what experimental techniques are applied, etc. Should be treated
    /// as a chunk of restructuredText (http://docutils.sourceforge.net/rst.html),
    /// with equations to be rendered as LaTeX with amsmath operators.
    virtual std::string description() const {
      return info().description();
    }

    /// @brief Information about the events needed as input for this analysis.
    ///
    /// Event types, energies, kinematic cuts, particles to be considered
    /// stable, etc. etc. Should be treated as a restructuredText bullet list
    /// (http://docutils.sourceforge.net/rst.html)
    virtual std::string runInfo() const {
      return info().runInfo();
    }

    /// Experiment which performed and published this analysis.
    virtual std::string experiment() const {
      return info().experiment();
    }

    /// Collider on which the experiment ran.
    virtual std::string collider() const {
      return info().collider();
    }

    /// When the original experimental analysis was published.
    virtual std::string year() const {
      return info().year();
    }

    /// The integrated luminosity in inverse femtobarn
    virtual double luminosityfb() const {
      return info().luminosityfb();
    }
    /// The integrated luminosity in inverse picobarn
    virtual double luminosity() const {
      return info().luminosity();
    }

    /// Journal, and preprint references.
    virtual std::vector<std::string> references() const {
      return info().references();
    }

    /// BibTeX citation key for this article.
    virtual std::string bibKey() const {
      return info().bibKey();
    }

    /// BibTeX citation entry for this article.
    virtual std::string bibTeX() const {
      return info().bibTeX();
    }

    /// Whether this analysis is trusted (in any way!)
    virtual std::string status() const {
      return (info().status().empty()) ? "UNVALIDATED" : info().status();
    }

    /// Any work to be done on this analysis.
    virtual std::vector<std::string> todos() const {
      return info().todos();
    }

    /// make-style commands for validating this analysis.
    virtual std::vector<std::string> validation() const {
      return info().validation();
    }

    /// Does this analysis have a reentrant finalize()?
    virtual bool reentrant() const {
      return info().reentrant();
    }


    /// Location of reference data YODA file
    virtual std::string refFile() const {
      return info().refFile();
    }


    /// Return the allowed pairs of incoming beams required by this analysis.
    virtual const std::vector<PdgIdPair>& requiredBeams() const {
      return info().beams();
    }
    /// Declare the allowed pairs of incoming beams required by this analysis.
    virtual Analysis& setRequiredBeams(const std::vector<PdgIdPair>& requiredBeams) {
      info().setBeams(requiredBeams);
      return *this;
    }

    /// Sets of valid beam energy pairs, in GeV
    virtual const std::vector<std::pair<double, double> >& requiredEnergies() const {
      return info().energies();
    }

    /// Get vector of analysis keywords
    virtual const std::vector<std::string> & keywords() const {
      return info().keywords();
    }

    /// Declare the list of valid beam energy pairs, in GeV
    virtual Analysis& setRequiredEnergies(const std::vector<std::pair<double, double> >& requiredEnergies) {
      info().setEnergies(requiredEnergies);
      return *this;
    }


    /// Get the actual AnalysisInfo object in which all this metadata is stored (non-const).
    /// @note For *internal* use!
    AnalysisInfo& info() {
      assert(_info && "No AnalysisInfo object :O");
      return *_info;
    }

    /// @}


    /// @defgroup analysis_run Run conditions
    /// @{

    /// Incoming beams for this run
    const ParticlePair& beams() const;

    /// Incoming beam IDs for this run
    const PdgIdPair beamIds() const;

    /// Centre of mass energy for this run
    double sqrtS() const;

    /// Check if we are running rivet-merge
    bool merging() const {
      return sqrtS() <= 0.0;
    }

    /// @}


    /// @defgroup analysis_beamcompat Analysis / beam compatibility testing
    ///
    /// @todo Replace with beamsCompatible() with no args (calling beams() function internally)
    /// @todo Add beamsMatch() methods with same (shared-code?) tolerance as in beamsCompatible()
    /// @{

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(const ParticlePair& beams) const;

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(PdgId beam1, PdgId beam2, double e1, double e2) const;

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(const PdgIdPair& beams, const std::pair<double,double>& energies) const;

    /// @}

    /// Access the controlling AnalysisHandler object.
    AnalysisHandler& handler() const { return *_analysishandler; }


  protected:

    /// Get a Log object based on the name() property of the calling analysis object.
    Log& getLog() const;

    /// Get the process cross-section in pb. Throws if this hasn't been set.
    double crossSection() const;

    /// Get the process cross-section per generated event in pb. Throws if this
    /// hasn't been set.
    double crossSectionPerEvent() const;

    /// Get the process cross-section error in pb. Throws if this hasn't been set.
    double crossSectionError() const;

    /// Get the process cross-section error per generated event in
    /// pb. Throws if this hasn't been set.
    double crossSectionErrorPerEvent() const;

    /// @brief Get the number of events seen (via the analysis handler).
    ///
    /// @note Use in the finalize phase only.
    size_t numEvents() const;

    /// @brief Get the sum of event weights seen (via the analysis handler).
    ///
    /// @note Use in the finalize phase only.
    double sumW() const;
    /// Alias
    double sumOfWeights() const { return sumW(); }

    /// @brief Get the sum of squared event weights seen (via the analysis handler).
    ///
    /// @note Use in the finalize phase only.
    double sumW2() const;


  protected:

    /// @defgroup analysis_histopaths Histogram paths
    ///
    /// @todo Add "tmp" flags to return /ANA/TMP/foo/bar paths
    ///
    /// @{

    /// Get the canonical histogram "directory" path for this analysis.
    const std::string histoDir() const;

    /// Get the canonical histogram path for the named histogram in this analysis.
    const std::string histoPath(const std::string& hname) const;

    /// Get the canonical histogram path for the numbered histogram in this analysis.
    const std::string histoPath(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const;

    /// Get the internal histogram name for given d, x and y (cf. HepData)
    const std::string mkAxisCode(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const;

    /// @}


    /// @defgroup analysis_refdata Histogram reference data
    /// @{

    /// Get all reference data objects for this analysis
    const std::map<std::string, YODA::AnalysisObjectPtr>& refData() const {
      _cacheRefData();
      return _refdata;
    }


    /// Get reference data for a named histo
    /// @todo SFINAE to ensure that the type inherits from YODA::AnalysisObject?
    template <typename T=YODA::Scatter2D>
    const T& refData(const string& hname) const {
      _cacheRefData();
      MSG_TRACE("Using histo bin edges for " << name() << ":" << hname);
      if (!_refdata[hname]) {
        MSG_ERROR("Can't find reference histogram " << hname);
        throw Exception("Reference data " + hname + " not found.");
      }
      return dynamic_cast<T&>(*_refdata[hname]);
    }


    /// Get reference data for a numbered histo
    /// @todo SFINAE to ensure that the type inherits from YODA::AnalysisObject?
    template <typename T=YODA::Scatter2D>
    const T& refData(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
      const string hname = mkAxisCode(datasetId, xAxisId, yAxisId);
      return refData<T>(hname);
    }

    /// @}


    /// @defgroup analysis_cbook Counter booking
    ///
    /// @todo Add "tmp" flags to book in standard temporary paths
    ///
    /// @{

    /// Book a counter.
    CounterPtr& book(CounterPtr&, const std::string& name);

    /// Book a counter, using a path generated from the dataset and axis ID codes
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    CounterPtr& book(CounterPtr&, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId);

    /// @}


    /// @defgroup analysis_h1book 1D histogram booking
    /// @{

    /// Book a 1D histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    Histo1DPtr& book(Histo1DPtr&,const std::string& name, size_t nbins, double lower, double upper);

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Histo1DPtr& book(Histo1DPtr&,const std::string& name, const std::vector<double>& binedges);

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Histo1DPtr& book(Histo1DPtr&,const std::string& name, const std::initializer_list<double>& binedges);

    /// Book a 1D histogram with binning from a reference scatter.
    Histo1DPtr& book(Histo1DPtr&,const std::string& name, const Scatter2D& refscatter);

    /// Book a 1D histogram, using the binnings in the reference data histogram.
    Histo1DPtr& book(Histo1DPtr&,const std::string& name);

    /// Book a 1D histogram, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    Histo1DPtr& book(Histo1DPtr&,unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId);

    /// @}


    /// @defgroup analysis_h2book 2D histogram booking
    /// @{

    /// Book a 2D histogram with @a nxbins and @a nybins uniformly
    /// distributed across the ranges @a xlower - @a xupper and @a
    /// ylower - @a yupper respectively along the x- and y-axis.
    Histo2DPtr& book(Histo2DPtr&,const std::string& name,
                           size_t nxbins, double xlower, double xupper,
                           size_t nybins, double ylower, double yupper);

    /// Book a 2D histogram with non-uniform bins defined by the
    /// vectors of bin edges @a xbinedges and @a ybinedges.
    Histo2DPtr& book(Histo2DPtr&,const std::string& name,
                           const std::vector<double>& xbinedges,
                           const std::vector<double>& ybinedges);

    /// Book a 2D histogram with non-uniform bins defined by the
    /// vectors of bin edges @a xbinedges and @a ybinedges.
    Histo2DPtr& book(Histo2DPtr&,const std::string& name,
                           const std::initializer_list<double>& xbinedges,
                           const std::initializer_list<double>& ybinedges);

    /// Book a 2D histogram with binning from a reference scatter.
    Histo2DPtr& book(Histo2DPtr&,const std::string& name,
                           const Scatter3D& refscatter);

    /// Book a 2D histogram, using the binnings in the reference data histogram.
    Histo2DPtr& book(Histo2DPtr&,const std::string& name);

    /// Book a 2D histogram, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    Histo2DPtr& book(Histo2DPtr&,unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId);

    /// @}


    /// @defgroup analysis_p1book 1D profile histogram booking
    /// @{

    /// Book a 1D profile histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    Profile1DPtr& book(Profile1DPtr&,  const std::string& name, size_t nbins, double lower, double upper);

    /// Book a 1D profile histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Profile1DPtr& book(Profile1DPtr&,  const std::string& name, const std::vector<double>& binedges);

    /// Book a 1D profile histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    Profile1DPtr& book(Profile1DPtr&,  const std::string& name, const std::initializer_list<double>& binedges);

    /// Book a 1D profile histogram with binning from a reference scatter.
    Profile1DPtr& book(Profile1DPtr&,  const std::string& name, const Scatter2D& refscatter);

    /// Book a 1D profile histogram, using the binnings in the reference data histogram.
    Profile1DPtr& book(Profile1DPtr&,  const std::string& name);

    /// Book a 1D profile histogram, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    Profile1DPtr& book(Profile1DPtr&,  unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId);

    /// @}


    /// @defgroup analysis_p2book 2D profile histogram booking
    /// @{

    /// Book a 2D profile histogram with @a nxbins and @a nybins uniformly
    /// distributed across the ranges @a xlower - @a xupper and @a ylower - @a
    /// yupper respectively along the x- and y-axis.
    Profile2DPtr& book(Profile2DPtr&,  const std::string& name,
                               size_t nxbins, double xlower, double xupper,
                               size_t nybins, double ylower, double yupper);

    /// Book a 2D profile histogram with non-uniform bins defined by the vectorx
    /// of bin edges @a xbinedges and @a ybinedges.
    Profile2DPtr& book(Profile2DPtr&,  const std::string& name,
                               const std::vector<double>& xbinedges,
                               const std::vector<double>& ybinedges);

    /// Book a 2D profile histogram with non-uniform bins defined by the vectorx
    /// of bin edges @a xbinedges and @a ybinedges.
    Profile2DPtr& book(Profile2DPtr&,  const std::string& name,
                               const std::initializer_list<double>& xbinedges,
                               const std::initializer_list<double>& ybinedges);

    /// @todo REINSTATE

    // /// Book a 2D profile histogram with binning from a reference scatter.
    // Profile2DPtr& book(const Profile2DPtr&, const std::string& name,
    //                            const Scatter3D& refscatter);

    // /// Book a 2D profile histogram, using the binnings in the reference data histogram.
    // Profile2DPtr& book(const Profile2DPtr&, const std::string& name);

    // /// Book a 2D profile histogram, using the binnings in the reference data histogram.
    // ///
    // /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    // Profile2DPtr& book(const Profile2DPtr&, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId);

    /// @}


    /// @defgroup analysis_s2book 2D scatter booking
    /// @{

    /// @brief Book a 2-dimensional data point set with the given name.
    ///
    /// @note Unlike histogram booking, scatter booking by default makes no
    /// attempt to use reference data to pre-fill the data object. If you want
    /// this, which is sometimes useful e.g. when the x-position is not really
    /// meaningful and can't be extracted from the data, then set the @a
    /// copy_pts parameter to true. This creates points to match the reference
    /// data's x values and errors, but with the y values and errors zeroed...
    /// assuming that there is a reference histo with the same name: if there
    /// isn't, an exception will be thrown.
    Scatter2DPtr& book(Scatter2DPtr& s2d, const string& hname, bool copy_pts = false);

    /// @brief Book a 2-dimensional data point set, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    ///
    /// @note Unlike histogram booking, scatter booking by default makes no
    /// attempt to use reference data to pre-fill the data object. If you want
    /// this, which is sometimes useful e.g. when the x-position is not really
    /// meaningful and can't be extracted from the data, then set the @a
    /// copy_pts parameter to true. This creates points to match the reference
    /// data's x values and errors, but with the y values and errors zeroed.
    Scatter2DPtr& book(Scatter2DPtr& s2d, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId, bool copy_pts = false);

    /// @brief Book a 2-dimensional data point set with equally spaced x-points in a range.
    ///
    /// The y values and errors will be set to 0.
    Scatter2DPtr& book(Scatter2DPtr& s2d, const string& hname, size_t npts, double lower, double upper);

    /// @brief Book a 2-dimensional data point set based on provided contiguous "bin edges".
    ///
    /// The y values and errors will be set to 0.
    Scatter2DPtr& book(Scatter2DPtr& s2d, const string& hname, const std::vector<double>& binedges);

    /// Book a 2-dimensional data point set with x-points from an existing scatter and a new path.
    Scatter2DPtr& book(Scatter2DPtr& s2d, const string& hname, const Scatter2D& refscatter);

    /// @}

    /// @defgroup analysis_s3book 3D scatter booking
    /// @{

    /// @brief Book a 3-dimensional data point set with the given name.
    ///
    /// @note Unlike histogram booking, scatter booking by default makes no
    /// attempt to use reference data to pre-fill the data object. If you want
    /// this, which is sometimes useful e.g. when the x-position is not really
    /// meaningful and can't be extracted from the data, then set the @a
    /// copy_pts parameter to true. This creates points to match the reference
    /// data's x values and errors, but with the y values and errors zeroed...
    /// assuming that there is a reference histo with the same name: if there
    /// isn't, an exception will be thrown.
    Scatter3DPtr& book(Scatter3DPtr& s3d, const std::string& hname, bool copy_pts=false);

    /// @brief Book a 3-dimensional data point set, using the binnings in the reference data histogram.
    ///
    /// The paper, dataset and x/y-axis IDs will be used to build the histo name in the HepData standard way.
    ///
    /// @note Unlike histogram booking, scatter booking by default makes no
    /// attempt to use reference data to pre-fill the data object. If you want
    /// this, which is sometimes useful e.g. when the x-position is not really
    /// meaningful and can't be extracted from the data, then set the @a
    /// copy_pts parameter to true. This creates points to match the reference
    /// data's x values and errors, but with the y values and errors zeroed.
    Scatter3DPtr& book(Scatter3DPtr& s3d, unsigned int datasetId, unsigned int xAxisId,
                        unsigned int yAxisId, unsigned int zAxisId, bool copy_pts=false);

    /// @brief Book a 3-dimensional data point set with equally spaced x-points in a range.
    ///
    /// The y values and errors will be set to 0.
    Scatter3DPtr& book(Scatter3DPtr& s3d, const std::string& hname,
                               size_t xnpts, double xlower, double xupper,
                               size_t ynpts, double ylower, double yupper);

    /// @brief Book a 3-dimensional data point set based on provided contiguous "bin edges".
    ///
    /// The y values and errors will be set to 0.
    Scatter3DPtr& book(Scatter3DPtr& s3d, const std::string& hname,
                               const std::vector<double>& xbinedges,
                               const std::vector<double>& ybinedges);

    /// Book a 3-dimensional data point set with x-points from an existing scatter and a new path.
    Scatter3DPtr& book(Scatter3DPtr& s3d, const std::string& hname, const Scatter3D& refscatter);

    /// @}


  public:

    /// @name Allow RAW histograms to be read in to local objects.
    /// @todo Should be protected, not public?
    /// @todo Why is the function body written this way? To avoid the virtual function being optimised away?
    virtual void rawHookIn(YODA::AnalysisObjectPtr yao) {
      (void) yao;
    }

    /// @name Provide access to RAW histograms before writing out to file.
    /// @todo Should be protected, not public?
    /// @todo Signature should pass the vector by reference?
    /// @todo Why is the function body written this way? To avoid the virtual function being optimised away?
    virtual void rawHookOut(vector<MultiweightAOPtr> raos, size_t iW) {
      (void) raos;
      (void) iW;
    }

    /// @name Accessing options for this Analysis instance.
    //@{

    /// Return the map of all options given to this analysis.
    const std::map<std::string,std::string>& options() const {
      return _options;
    }

    /// Get an option for this analysis instance as a string.
    std::string getOption(std::string optname) const {
      if ( _options.find(optname) != _options.end() )
        return _options.find(optname)->second;
      return "";
    }

    /// @brief Get an option for this analysis instance converted to a specific type
    ///
    /// The return type is given by the specified @a def value, or by an explicit template
    /// type-argument, e.g. getOption<double>("FOO", 3).
    template<typename T>
    T getOption(std::string optname, T def) const {
      if (_options.find(optname) == _options.end()) return def;
      std::stringstream ss;
      ss << _options.find(optname)->second;
      T ret;
      ss >> ret;
      return ret;
    }

    /// @brief Sane overload for literal character strings (which don't play well with stringstream)
    ///
    /// Note this isn't a template specialisation, because we can't return a non-static
    /// char*, and T-as-return-type is built into the template function definition.
    std::string getOption(std::string optname, const char* def) {
      return getOption<std::string>(optname, def);
    }

    /// @}


    /// @defgroup analysis_bookhi Booking heavy ion features
    /// @{

    /// @brief Book a CentralityProjection
    ///
    /// Using a SingleValueProjection, @a proj, giving the value of an
    /// experimental observable to be used as a centrality estimator,
    /// book a CentralityProjection based on the experimentally
    /// measured pecentiles of this observable (as given by the
    /// reference data for the @a calHistName histogram in the @a
    /// calAnaName analysis. If a preloaded file with the output of a
    /// run using the @a calAnaName analysis contains a valid
    /// generated @a calHistName histogram, it will be used as an
    /// optional percentile binning. Also if this preloaded file
    /// contains a histogram with the name @a calHistName with an
    /// appended "_IMP" This histogram will be used to add an optional
    /// centrality percentile based on the generated impact
    /// parameter. If @a increasing is true, a low (high) value of @a proj
    /// is assumed to correspond to a more peripheral (central) event.
    const CentralityProjection&
    declareCentrality(const SingleValueProjection &proj,
                      string calAnaName, string calHistName,
                      const string projName, bool increasing=false);


    /// @brief Book a Percentile wrapper around AnalysisObjects.
    ///
    /// Based on a previously registered CentralityProjection named @a
    /// projName book one AnalysisObject for each @a centralityBin and
    /// name them according to the corresponding code in the @a ref
    /// vector.
    ///
    /// @todo Convert to just be called book() cf. others
    template <class T>
    Percentile<T> bookPercentile(string projName,
                                 vector<pair<float, float> > centralityBins,
                                 vector<tuple<int, int, int> > ref) {

      typedef typename ReferenceTraits<T>::RefT RefT;
      typedef rivet_shared_ptr<Wrapper<T>> WrapT;

      Percentile<T> pctl(this, projName);

      const int nCent = centralityBins.size();
      for (int iCent = 0; iCent < nCent; ++iCent) {
        const string axisCode = mkAxisCode(std::get<0>(ref[iCent]),
                                           std::get<1>(ref[iCent]),
                                           std::get<2>(ref[iCent]));
        const RefT & refscatter = refData<RefT>(axisCode);

        WrapT wtf(_weightNames(), T(refscatter, histoPath(axisCode)));
        wtf = addAnalysisObject(wtf);

        CounterPtr cnt(_weightNames(), Counter(histoPath("TMP/COUNTER/" + axisCode)));
        cnt = addAnalysisObject(cnt);

        pctl.add(wtf, cnt, centralityBins[iCent]);
      }
      return pctl;
    }


    // /// @brief Book Percentile wrappers around AnalysisObjects.
    // ///
    // /// Based on a previously registered CentralityProjection named @a
    // /// projName book one (or several) AnalysisObject(s) named
    // /// according to @a ref where the x-axis will be filled according
    // /// to the percentile output(s) of the @projName.
    // ///
    // /// @todo Convert to just be called book() cf. others
    // template <class T>
    // PercentileXaxis<T> bookPercentileXaxis(string projName,
    //                                        tuple<int, int, int> ref) {

    //   typedef typename ReferenceTraits<T>::RefT RefT;
    //   typedef rivet_shared_ptr<Wrapper<T>> WrapT;

    //   PercentileXaxis<T> pctl(this, projName);

    //   const string axisCode = mkAxisCode(std::get<0>(ref),
    //                                      std::get<1>(ref),
    //                                      std::get<2>(ref));
    //   const RefT & refscatter = refData<RefT>(axisCode);

    //   WrapT wtf(_weightNames(), T(refscatter, histoPath(axisCode)));
    //   wtf = addAnalysisObject(wtf);

    //   CounterPtr cnt(_weightNames(), Counter());
    //   cnt = addAnalysisObject(cnt);

    //   pctl.add(wtf, cnt);
    //   return pctl;
    // }

    /// @}


  private:

    // Functions that have to be defined in the .cc file to avoid circular #includes

    /// Get the list of weight names from the handler
    vector<string> _weightNames() const;

    /// Get the list of weight names from the handler
    YODA::AnalysisObjectPtr _getPreload(string name) const;

    /// Get an AO from another analysis
    MultiweightAOPtr _getOtherAnalysisObject(const std::string & ananame, const std::string& name);

    /// Check that analysis objects aren't being booked/registered outside the init stage
    void _checkBookInit() const;

    /// Check if we are in the init stage.
    bool _inInit() const;

    /// Check if we are in the finalize stage.
    bool _inFinalize() const;


  private:

    /// To be used in finalize context only
    class CounterAdapter {
    public:

      CounterAdapter(double x) : x_(x) {}

      CounterAdapter(const YODA::Counter & c) : x_(c.val()) {}

      CounterAdapter(const YODA::Scatter1D & s) : x_(s.points()[0].x()) {
        assert( s.numPoints() == 1 || "Can only scale by a single value.");
      }

      operator double() const { return x_; }

    private:
      double x_;

    };


  public:

    double dbl(double          x) { return x; }
    double dbl(const YODA::Counter   & c) { return c.val(); }
    double dbl(const YODA::Scatter1D & s) {
      assert( s.numPoints() == 1 );
      return s.points()[0].x();
    }


    /// @defgroup analysis_manip Analysis object manipulation
    ///
    /// @todo Should really be protected: only public to keep BinnedHistogram happy for now...
    /// @{

    /// @todo Add functions for converting histos and profiles to scatters (with bin-width and focus control, over syst variations)

    /// Multiplicatively scale the given counter, @a cnt, by factor @a factor.
    void scale(CounterPtr cnt, CounterAdapter factor);

    /// Multiplicatively scale the given counters, @a cnts, by factor @a factor.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of CounterPtrs
    void scale(const std::vector<CounterPtr>& cnts, CounterAdapter factor) {
      for (auto& c : cnts) scale(c, factor);
    }

    /// Iteratively scale the counters in the map @a maps, by factor @a factor.
    template<typename T>
    void scale(const std::map<T, CounterPtr>& maps, CounterAdapter factor) {
      for (auto& m : maps) scale(m.second, factor);
    }

    /// @todo YUCK!
    template <std::size_t array_size>
    void scale(const CounterPtr (&cnts)[array_size], CounterAdapter factor) {
      // for (size_t i = 0; i < std::extent<decltype(cnts)>::value; ++i) scale(cnts[i], factor);
      for (auto& c : cnts) scale(c, factor);
    }


    /// Normalize the given histogram, @a histo, to area = @a norm.
    void normalize(Histo1DPtr histo, CounterAdapter norm=1.0, bool includeoverflows=true);

    /// Normalize the given histograms, @a histos, to area = @a norm.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo1DPtrs
    void normalize(const std::vector<Histo1DPtr>& histos, CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }

    /// Normalize the histograms in map, @a maps, to area = @a norm.
    template<typename T>
    void normalize(const std::map<T, Histo1DPtr>& maps, CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& m : maps) normalize(m.second, norm, includeoverflows);
    }

    /// @todo YUCK!
    template <std::size_t array_size>
    void normalize(const Histo1DPtr (&histos)[array_size], CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }

    /// Multiplicatively scale the given histogram, @a histo, by factor @a factor.
    void scale(Histo1DPtr histo, CounterAdapter factor);

    /// Multiplicatively scale the given histograms, @a histos, by factor @a factor.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo1DPtrs
    void scale(const std::vector<Histo1DPtr>& histos, CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }

    /// Iteratively scale the histograms in the map, @a maps, by factor @a factor.
    template<typename T>
    void scale(const std::map<T, Histo1DPtr>& maps, CounterAdapter factor) {
      for (auto& m : maps) scale(m.second, factor);
    }

    /// @todo YUCK!
    template <std::size_t array_size>
    void scale(const Histo1DPtr (&histos)[array_size], CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }


    /// Normalize the given histogram, @a histo, to area = @a norm.
    void normalize(Histo2DPtr histo, CounterAdapter norm=1.0, bool includeoverflows=true);

    /// Normalize the given histograms, @a histos, to area = @a norm.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo2DPtrs
    void normalize(const std::vector<Histo2DPtr>& histos, CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }

    /// Normalize the histograms in map, @a maps, to area = @a norm.
    template<typename T>
    void normalize(const std::map<T, Histo2DPtr>& maps, CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& m : maps) normalize(m.second, norm, includeoverflows);
    }

    /// @todo YUCK!
    template <std::size_t array_size>
    void normalize(const Histo2DPtr (&histos)[array_size], CounterAdapter norm=1.0, bool includeoverflows=true) {
      for (auto& h : histos) normalize(h, norm, includeoverflows);
    }

    /// Multiplicatively scale the given histogram, @a histo, by factor @a factor.
    void scale(Histo2DPtr histo, CounterAdapter factor);

    /// Multiplicatively scale the given histograms, @a histos, by factor @a factor.
    /// @note Constness intentional, if weird, to allow passing rvalue refs of smart ptrs (argh)
    /// @todo Use SFINAE for a generic iterable of Histo2DPtrs
    void scale(const std::vector<Histo2DPtr>& histos, CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }

    /// Iteratively scale the histograms in the map, @a maps, by factor @a factor.
    template<typename T>
    void scale(const std::map<T, Histo2DPtr>& maps, CounterAdapter factor) {
      for (auto& m : maps) scale(m.second, factor);
    }

    /// @todo YUCK!
    template <std::size_t array_size>
    void scale(const Histo2DPtr (&histos)[array_size], CounterAdapter factor) {
      for (auto& h : histos) scale(h, factor);
    }


    /// Helper for counter division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(CounterPtr c1, CounterPtr c2, Scatter1DPtr s) const;

    /// Helper for histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Counter& c1, const YODA::Counter& c2, Scatter1DPtr s) const;


    /// Helper for histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const;

    /// Helper for histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Histo1D& h1, const YODA::Histo1D& h2, Scatter2DPtr s) const;


    /// Helper for profile histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Profile1DPtr p1, Profile1DPtr p2, Scatter2DPtr s) const;

    /// Helper for profile histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Profile1D& p1, const YODA::Profile1D& p2, Scatter2DPtr s) const;


    /// Helper for 2D histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Histo2DPtr h1, Histo2DPtr h2, Scatter3DPtr s) const;

    /// Helper for 2D histogram division with raw YODA objects.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(const YODA::Histo2D& h1, const YODA::Histo2D& h2, Scatter3DPtr s) const;


    /// Helper for 2D profile histogram division.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void divide(Profile2DPtr p1, Profile2DPtr p2, Scatter3DPtr s) const;

    /// Helper for 2D profile histogram division with raw YODA objects
    ///
    /// @note Assigns to the (already registered) output scatter, @a s.  Preserves the path information of the target.
    void divide(const YODA::Profile2D& p1, const YODA::Profile2D& p2, Scatter3DPtr s) const;


    /// Helper for histogram efficiency calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void efficiency(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const;

    /// Helper for histogram efficiency calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void efficiency(const YODA::Histo1D& h1, const YODA::Histo1D& h2, Scatter2DPtr s) const;


    /// Helper for histogram asymmetry calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void asymm(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const;

    /// Helper for histogram asymmetry calculation.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void asymm(const YODA::Histo1D& h1, const YODA::Histo1D& h2, Scatter2DPtr s) const;


    /// Helper for converting a differential histo to an integral one.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void integrate(Histo1DPtr h, Scatter2DPtr s) const;

    /// Helper for converting a differential histo to an integral one.
    ///
    /// @note Assigns to the (already registered) output scatter, @a s. Preserves the path information of the target.
    void integrate(const Histo1D& h, Scatter2DPtr s) const;

    /// @}


  public:

    /// List of registered analysis data objects
    const vector<MultiweightAOPtr>& analysisObjects() const {
      return _analysisobjects;
    }


  protected:

    /// @defgroup analysis_aoaccess Data object registration, retrieval, and removal
    /// @{

    /// Get the default/nominal weight index
    size_t defaultWeightIndex() const;

    /// Get a preloaded YODA object.
    template <typename YODAT>
    shared_ptr<YODAT> getPreload(string path) const {
      return dynamic_pointer_cast<YODAT>(_getPreload(path));
    }


    /// Register a new data object, optionally read in preloaded data.
    template <typename YODAT>
    rivet_shared_ptr< Wrapper<YODAT> > registerAO(const YODAT& yao) {
      typedef Wrapper<YODAT> WrapperT;
      typedef shared_ptr<YODAT> YODAPtrT;
      typedef rivet_shared_ptr<WrapperT> RAOT;

      if ( !_inInit() && !_inFinalize() ) {
        MSG_ERROR("Can't book objects outside of init() or finalize()");
        throw UserError(name() + ": Can't book objects outside of init() or finalize().");
      }

      // First check that we haven't booked this before.
      // This is allowed when booking in finalize: just warn in that case.
      // If in init(), throw an exception: it's 99.9% never going to be intentional.
      for (auto& waold : analysisObjects()) {
        if ( yao.path() == waold.get()->basePath() ) {
          const string msg = "Found double-booking of " + yao.path() + " in " + name();
          if ( _inInit() ) {
            MSG_ERROR(msg);
            throw LookupError(msg);
          } else {
            MSG_WARNING(msg + ". Keeping previous booking");
          }
          return RAOT(dynamic_pointer_cast<WrapperT>(waold.get()));
        }
      }

      shared_ptr<WrapperT> wao = make_shared<WrapperT>();
      wao->_basePath = yao.path();
      YODAPtrT yaop = make_shared<YODAT>(yao);

      for (const string& weightname : _weightNames()) {
        // Create two YODA objects for each weight. Copy from
        // preloaded YODAs if present. First the finalized yoda:
        string finalpath = yao.path();
        if ( weightname != "" ) finalpath +=  "[" + weightname + "]";
        YODAPtrT preload = getPreload<YODAT>(finalpath);
        if ( preload ) {
          if ( !bookingCompatible(preload, yaop) ) {
            /// @todo What about if/when we want to make the final objects the Scatter or binned persistent type?
            MSG_WARNING("Found incompatible pre-existing data object with same base path "
                        << finalpath <<  " for " << name());
            preload = nullptr;
          } else {
            MSG_TRACE("Using preloaded " << finalpath << " in " <<name());
            wao->_final.push_back(make_shared<YODAT>(*preload));
          }
        }
        if ( !preload ) {
          wao->_final.push_back(make_shared<YODAT>(yao));
          wao->_final.back()->setPath(finalpath);
        }

        // Then the raw filling yodas.
        string rawpath = "/RAW" + finalpath;
        preload = getPreload<YODAT>(rawpath);
        if ( preload ) {
          if ( !bookingCompatible(preload, yaop) ) {
            MSG_WARNING("Found incompatible pre-existing data object with same base path "
                        << rawpath <<  " for " << name());
            preload = nullptr;
          } else {
            MSG_TRACE("Using preloaded " << rawpath << " in " <<name());
            wao->_persistent.push_back(make_shared<YODAT>(*preload));
          }
        }
        if ( !preload ) {
          wao->_persistent.push_back(make_shared<YODAT>(yao));
          wao->_persistent.back()->setPath(rawpath);
        }
      }
      rivet_shared_ptr<WrapperT> ret(wao);

      ret.get()->unsetActiveWeight();
      if ( _inFinalize() ) {
        // If booked in finalize() we assume it is the first time
        // finalize is run.
        ret.get()->pushToFinal();
        ret.get()->setActiveFinalWeightIdx(0);
      }
      _analysisobjects.push_back(ret);

      return ret;
    }


    /// Register a data object in the histogram system
    template <typename AO=MultiweightAOPtr>
    AO addAnalysisObject(const AO& aonew) {
      _checkBookInit();

      for (const MultiweightAOPtr& ao : analysisObjects()) {

        // Check AO base-name first
        ao.get()->setActiveWeightIdx(defaultWeightIndex());
        aonew.get()->setActiveWeightIdx(defaultWeightIndex());
        if (ao->path() != aonew->path()) continue;

        // If base-name matches, check compatibility
        // NB. This evil is because dynamic_ptr_cast can't work on rivet_shared_ptr directly
        AO aoold = AO(dynamic_pointer_cast<typename AO::value_type>(ao.get())); //< OMG
        if ( !aoold || !bookingCompatible(aonew, aoold) ) {
          MSG_WARNING("Found incompatible pre-existing data object with same base path "
                      << aonew->path() <<  " for " << name());
          throw LookupError("Found incompatible pre-existing data object with same base path during AO booking");
        }

        // Finally, check all weight variations
        for (size_t weightIdx = 0; weightIdx < _weightNames().size(); ++weightIdx) {
          aoold.get()->setActiveWeightIdx(weightIdx);
          aonew.get()->setActiveWeightIdx(weightIdx);
          if (aoold->path() != aonew->path()) {
            MSG_WARNING("Found incompatible pre-existing data object with different weight-path "
                        << aonew->path() <<  " for " << name());
            throw LookupError("Found incompatible pre-existing data object with same weight-path during AO booking");
          }
        }

        // They're fully compatible: bind and return
        aoold.get()->unsetActiveWeight();
        MSG_TRACE("Bound pre-existing data object " << aoold->path() <<  " for " << name());
        return aoold;
      }

      // No equivalent found
      MSG_TRACE("Registered " << aonew->annotation("Type") << " " << aonew->path() <<  " for " << name());
      aonew.get()->unsetActiveWeight();

      _analysisobjects.push_back(aonew);
      return aonew;
    }

    /// Unregister a data object from the histogram system (by name)
    void removeAnalysisObject(const std::string& path);

    /// Unregister a data object from the histogram system (by pointer)
    void removeAnalysisObject(const MultiweightAOPtr& ao);

    // /// Get all data objects, for all analyses, from the AnalysisHandler
    // /// @todo Can we remove this? Why not call handler().getData()?
    // vector<YODA::AnalysisObjectPtr> getAllData(bool includeorphans) const;


    /// Get a Rivet data object from the histogram system
    template <typename AO=MultiweightAOPtr>
    const AO getAnalysisObject(const std::string& aoname) const {
      for (const MultiweightAOPtr& ao : analysisObjects()) {
        ao.get()->setActiveWeightIdx(defaultWeightIndex());
        if (ao->path() == histoPath(aoname)) {
          // return dynamic_pointer_cast<AO>(ao);
          return AO(dynamic_pointer_cast<typename AO::value_type>(ao.get()));
        }
      }
      throw LookupError("Data object " + histoPath(aoname) + " not found");
    }


    // /// Get a data object from the histogram system
    // template <typename AO=YODA::AnalysisObject>
    // const std::shared_ptr<AO> getAnalysisObject(const std::string& name) const {
    //   foreach (const AnalysisObjectPtr& ao, analysisObjects()) {
    //     if (ao->path() == histoPath(name)) return dynamic_pointer_cast<AO>(ao);
    //   }
    //   throw LookupError("Data object " + histoPath(name) + " not found");
    // }

    // /// Get a data object from the histogram system (non-const)
    // template <typename AO=YODA::AnalysisObject>
    // std::shared_ptr<AO> getAnalysisObject(const std::string& name) {
    //   foreach (const AnalysisObjectPtr& ao, analysisObjects()) {
    //     if (ao->path() == histoPath(name)) return dynamic_pointer_cast<AO>(ao);
    //   }
    //   throw LookupError("Data object " + histoPath(name) + " not found");
    // }


    /// Get a data object from another analysis (e.g. preloaded
    /// calibration histogram).
    template <typename AO=MultiweightAOPtr>
    AO getAnalysisObject(const std::string& ananame,
                         const std::string& aoname) {
      MultiweightAOPtr ao = _getOtherAnalysisObject(ananame, aoname);
      // return dynamic_pointer_cast<AO>(ao);
      return AO(dynamic_pointer_cast<typename AO::value_type>(ao.get()));
    }


    // /// Get a named Histo1D object from the histogram system
    // const Histo1DPtr getHisto1D(const std::string& name) const {
    //   return getAnalysisObject<Histo1D>(name);
    // }

    // /// Get a named Histo1D object from the histogram system (non-const)
    // Histo1DPtr getHisto1D(const std::string& name) {
    //   return getAnalysisObject<Histo1D>(name);
    // }

    // /// Get a Histo1D object from the histogram system by axis ID codes (non-const)
    // const Histo1DPtr getHisto1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    //   return getAnalysisObject<Histo1D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }

    // /// Get a Histo1D object from the histogram system by axis ID codes (non-const)
    // Histo1DPtr getHisto1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    //   return getAnalysisObject<Histo1D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }


    // /// Get a named Histo2D object from the histogram system
    // const Histo2DPtr getHisto2D(const std::string& name) const {
    //   return getAnalysisObject<Histo2D>(name);
    // }

    // /// Get a named Histo2D object from the histogram system (non-const)
    // Histo2DPtr getHisto2D(const std::string& name) {
    //   return getAnalysisObject<Histo2D>(name);
    // }

    // /// Get a Histo2D object from the histogram system by axis ID codes (non-const)
    // const Histo2DPtr getHisto2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    //   return getAnalysisObject<Histo2D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }

    // /// Get a Histo2D object from the histogram system by axis ID codes (non-const)
    // Histo2DPtr getHisto2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    //   return getAnalysisObject<Histo2D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }


    // /// Get a named Profile1D object from the histogram system
    // const Profile1DPtr getProfile1D(const std::string& name) const {
    //   return getAnalysisObject<Profile1D>(name);
    // }

    // /// Get a named Profile1D object from the histogram system (non-const)
    // Profile1DPtr getProfile1D(const std::string& name) {
    //   return getAnalysisObject<Profile1D>(name);
    // }

    // /// Get a Profile1D object from the histogram system by axis ID codes (non-const)
    // const Profile1DPtr getProfile1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    //   return getAnalysisObject<Profile1D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }

    // /// Get a Profile1D object from the histogram system by axis ID codes (non-const)
    // Profile1DPtr getProfile1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    //   return getAnalysisObject<Profile1D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }


    // /// Get a named Profile2D object from the histogram system
    // const Profile2DPtr getProfile2D(const std::string& name) const {
    //   return getAnalysisObject<Profile2D>(name);
    // }

    // /// Get a named Profile2D object from the histogram system (non-const)
    // Profile2DPtr getProfile2D(const std::string& name) {
    //   return getAnalysisObject<Profile2D>(name);
    // }

    // /// Get a Profile2D object from the histogram system by axis ID codes (non-const)
    // const Profile2DPtr getProfile2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    //   return getAnalysisObject<Profile2D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }

    // /// Get a Profile2D object from the histogram system by axis ID codes (non-const)
    // Profile2DPtr getProfile2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    //   return getAnalysisObject<Profile2D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }


    // /// Get a named Scatter2D object from the histogram system
    // const Scatter2DPtr getScatter2D(const std::string& name) const {
    //   return getAnalysisObject<Scatter2D>(name);
    // }

    // /// Get a named Scatter2D object from the histogram system (non-const)
    // Scatter2DPtr getScatter2D(const std::string& name) {
    //   return getAnalysisObject<Scatter2D>(name);
    // }

    // /// Get a Scatter2D object from the histogram system by axis ID codes (non-const)
    // const Scatter2DPtr getScatter2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    //   return getAnalysisObject<Scatter2D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }

    // /// Get a Scatter2D object from the histogram system by axis ID codes (non-const)
    // Scatter2DPtr getScatter2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    //   return getAnalysisObject<Scatter2D>(makeAxisCode(datasetId, xAxisId, yAxisId));
    // }

    /// @}


  private:

    /// Name passed to constructor (used to find .info analysis data file, and as a fallback)
    string _defaultname;

    /// Pointer to analysis metadata object
    unique_ptr<AnalysisInfo> _info;

    /// Storage of all plot objects
    /// @todo Make this a map for fast lookup by path?
    vector<MultiweightAOPtr> _analysisobjects;

    /// @defgroup analysis_xsecvars Cross-section variables
    /// @{
    double _crossSection;
    bool _gotCrossSection;
    /// @}

    /// The controlling AnalysisHandler object.
    AnalysisHandler* _analysishandler;

    /// Collection of cached refdata to speed up many autobookings: the
    /// reference data file should only be read once.
    mutable std::map<std::string, YODA::AnalysisObjectPtr> _refdata;

    /// Options the (this instance of) the analysis
    map<string, string> _options;

    /// The string of options.
    string _optstring;


  private:

    /// @name Utility functions
    /// @{

    /// Get the reference data for this paper and cache it.
    void _cacheRefData() const;

    /// @}

  };


  // // Template specialisation for literal character strings (which don't play well with stringstream)
  // template<>
  // inline std::string Analysis::getOption(std::string optname, const char* def) {
  //   return getOption<std::string>(optname, def); //.c_str();
  // }


}


// Include definition of analysis plugin system so that analyses automatically see it when including Analysis.hh
#include "Rivet/AnalysisBuilder.hh"


/// @defgroup anamacros Analysis macros
/// @{

/// @def RIVET_DECLARE_PLUGIN
/// Preprocessor define to prettify the global-object plugin hook mechanism
#define RIVET_DECLARE_PLUGIN(clsname) ::Rivet::AnalysisBuilder<clsname> plugin_ ## clsname

/// @def RIVET_DECLARE_ALIASED_PLUGIN
/// Preprocessor define to prettify the global-object plugin hook mechanism, with an extra alias name for this analysis
#define RIVET_DECLARE_ALIASED_PLUGIN(clsname, alias) RIVET_DECLARE_PLUGIN(clsname)( #alias )

/// @def RIVET_DEFAULT_ANALYSIS_CTOR
/// Preprocessor define to prettify the awkward constructor with name string argument
#define RIVET_DEFAULT_ANALYSIS_CTOR(clsname) clsname() : Analysis(# clsname) {}



/// @def DECLARE_RIVET_PLUGIN
/// Preprocessor define to prettify the global-object plugin hook mechanism
///
/// @deprecated Prefer the RIVET_DECLARE_PLUGIN version with predictable RIVET_ prefix
#define DECLARE_RIVET_PLUGIN(clsname) ::Rivet::AnalysisBuilder<clsname> plugin_ ## clsname

/// @def DECLARE_ALIASED_RIVET_PLUGIN
/// Preprocessor define to prettify the global-object plugin hook mechanism, with an extra alias name for this analysis
///
/// @deprecated Prefer the RIVET_DECLARE_ALIASED_PLUGIN version with predictable RIVET_ prefix
// #define DECLARE_ALIASED_RIVET_PLUGIN(clsname, alias) Rivet::AnalysisBuilder<clsname> plugin_ ## clsname ## ( ## #alias ## )
#define DECLARE_ALIASED_RIVET_PLUGIN(clsname, alias) DECLARE_RIVET_PLUGIN(clsname)( #alias )

/// @def DEFAULT_RIVET_ANALYSIS_CONSTRUCTOR
/// Preprocessor define to prettify the awkward constructor with name string argument
///
/// @deprecated Prefer the "CTOR" version
#define DEFAULT_RIVET_ANALYSIS_CONSTRUCTOR(clsname) clsname() : Analysis(# clsname) {}

/// @def DEFAULT_RIVET_ANALYSIS_CTOR
/// Preprocessor define to prettify the awkward constructor with name string argument
///
/// @deprecated Prefer the RIVET_DEFAULT_ANALYSIS_CTOR version with predictable RIVET_ prefix
#define DEFAULT_RIVET_ANALYSIS_CTOR(clsname) DEFAULT_RIVET_ANALYSIS_CONSTRUCTOR(clsname)

/// @}


#endif

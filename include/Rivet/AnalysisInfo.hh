// -*- C++ -*-
#ifndef RIVET_AnalysisInfo_HH
#define RIVET_AnalysisInfo_HH

#include "Rivet/Config/RivetCommon.hh"
#include <ostream>

namespace Rivet {


  /// @brief Holder of analysis metadata
  class AnalysisInfo {
  public:

    /// Static factory method: returns null pointer if no metadata found
    static unique_ptr<AnalysisInfo> make(const std::string& name);

    /// Default constructor
    AnalysisInfo() { clear(); }

    ~AnalysisInfo() { }


    /// @defgroup anainfo_metadata Metadata
    ///
    /// Metadata is used for querying from the command line and also for
    /// building web pages and the analysis pages in the Rivet manual.
    /// @{

    /// Get the name of the analysis. By default this is computed using the
    /// experiment, year and Inspire/Spires ID metadata methods.
    std::string name() const {
      if (!_name.empty()) return _name;
      if (!experiment().empty() && !year().empty()) {
        if (!inspireId().empty()) {
          return experiment() + "_" + year() + "_I" + inspireId();
        } else if (!spiresId().empty()) {
          return experiment() + "_" + year() + "_S" + spiresId();
        }
      }
      return "";
    }
    /// Set the name of the analysis.
    void setName(const std::string& name) { _name = name; }

    /// Get the reference data name of the analysis (if different from plugin name).
    std::string getRefDataName() const {
      if (!_refDataName.empty())  return _refDataName;
      return name();
    }
    /// Set the reference data name of the analysis (if different from plugin name).
    void setRefDataName(const std::string& name) { _refDataName = name; }

    /// Get the Inspire (SPIRES replacement) ID code for this analysis.
    const std::string& inspireId() const { return _inspireId; }
    /// Set the Inspire (SPIRES replacement) ID code for this analysis.
    void setInspireId(const std::string& inspireId) { _inspireId = inspireId; }

    /// Get the SPIRES ID code for this analysis.
    ///
    /// @deprecated SPIRES itself is gone. Remove when all analyses migrated to I-names.
    const std::string& spiresId() const { return _spiresId; }
    /// Set the SPIRES ID code for this analysis.
    ///
    /// @deprecated SPIRES itself is gone. Remove when all analyses migrated to I-names.
    void setSpiresId(const std::string& spiresId) { _spiresId = spiresId; }

    /// @brief Names & emails of paper/analysis authors.
    ///
    /// Names and email of authors in 'NAME \<EMAIL\>' format. The first
    /// name in the list should be the primary contact person.
    const std::vector<std::string>& authors() const { return _authors; }
    /// Set the author list.
    void setAuthors(const std::vector<std::string>& authors) { _authors = authors; }

    /// @brief Get a short description of the analysis.
    ///
    /// Short (one sentence) description used as an index entry.
    /// Use @a description() to provide full descriptive paragraphs
    /// of analysis details.
    const std::string& summary() const { return _summary; }
    /// Set the short description for this analysis.
    void setSummary(const std::string& summary) { _summary = summary; }

    /// @brief Get a full description of the analysis.
    ///
    /// Full textual description of this analysis, what it is useful for,
    /// what experimental techniques are applied, etc. Should be treated
    /// as a chunk of restructuredText (http://docutils.sourceforge.net/rst.html),
    /// with equations to be rendered as LaTeX with amsmath operators.
    const std::string& description() const { return _description; }
    /// Set the full description for this analysis.
    void setDescription(const std::string& description) { _description = description; }

    /// @brief Information about the events needed as input for this analysis.
    ///
    /// Event types, energies, kinematic cuts, particles to be considered
    /// stable, etc. etc. Should be treated as a restructuredText bullet list
    /// (http://docutils.sourceforge.net/rst.html)
    const std::string& runInfo() const { return _runInfo; }
    /// Set the full description for this analysis.
    void setRunInfo(const std::string& runInfo) { _runInfo = runInfo; }

    /// Beam particle types
    const std::vector<PdgIdPair>& beams() const { return _beams; }
    /// Set beam particle types
    void setBeams(const std::vector<PdgIdPair>& beams) { _beams = beams; }

    /// Sets of valid beam energies
    const std::vector<std::pair<double,double> >& energies() const { return _energies; }
    /// Set the valid beam energies
    void setEnergies(const std::vector<std::pair<double, double> >& energies) { _energies = energies; }

    /// Experiment which performed and published this analysis.
    const std::string& experiment() const { return _experiment; }
    /// Set the experiment which performed and published this analysis.
    void setExperiment(const std::string& experiment) { _experiment = experiment; }

    /// Collider on which the experiment ran.
    const std::string& collider() const { return _collider; }
    /// Set the collider on which the experiment ran.
    void setCollider(const std::string& collider) { _collider = collider; }

    /// @brief When the original experimental analysis was published.
    ///
    /// When the refereed paper on which this is based was published, according to Inspire-HEP.
    const std::string& year() const { return _year; }

    /// Set the year in which the original experimental analysis was published.
    void setYear(const std::string& year) { _year = year; }

    /// The integrated data luminosity of the data set in 1/fb
    double luminosityfb() const { return _luminosityfb; }
    /// The integrated data luminosity of the data set in 1/pb
    double luminosity() const { return 1000*_luminosityfb; }

    /// Set the integrated data luminosity of the data set
    void setLuminosityfb(const double luminosityfb) { _luminosityfb = luminosityfb; }

    /// Journal and preprint references.
    const std::vector<std::string>& references() const { return _references; }
    /// Set the journal and preprint reference list.
    void setReferences(const std::vector<std::string>& references) { _references = references; }

    /// Analysis keywords, for grouping etc.
    const std::vector<std::string>& keywords() const { return _keywords; }
    void setKeywords(const std::vector<std::string>& keywords) { _keywords = keywords; }

    /// Any warning message
    const std::string& warning() const { return _warning; }
    void setWarning(const std::string warning) { _warning = warning; }

    /// Positive filtering regex for ref-data HepData sync
    const std::string& refMatch() const { return _refmatch; }
    void setRefMatch(const std::string refmatch) { _refmatch = refmatch; }

    /// Negative filtering regex for ref-data HepData sync
    const std::string& refUnmatch() const { return _refunmatch; }
    void setRefUnmatch(const std::string refunmatch) { _refunmatch = refunmatch; }

    /// Positive filtering regex for setting double precision in Writer
    const std::string& writerDoublePrecision() const { return _writerdoubleprecision; }
    void setWriterDoublePrecision(const std::string dp) { _writerdoubleprecision = dp; }

    /// BibTeX citation key for this article.
    const std::string& bibKey() const { return _bibKey; }
    /// Set the BibTeX citation key for this article.
    void setBibKey(const std::string& bibKey) { _bibKey = bibKey; }

    /// BibTeX citation entry for this article.
    const std::string& bibTeX() const { return _bibTeX; }
    /// Set the BibTeX citation entry for this article.
    void setBibTeX(const std::string& bibTeX) { _bibTeX = bibTeX; }

    /// Any work to be done on this analysis.
    const std::vector<std::string>& todos() const { return _todos; }
    /// Set the to-do list.
    void setTodos(const std::vector<std::string>& todos) { _todos = todos; }

    /// @}


    /// @defgroup anainfo_options Analysis-options support
    /// @{

    /// Get the option list.
    const std::vector<std::string>& options() const { return _options; }

    /// Check if the given option is valid.
    bool validOption(std::string key, std::string val) const;

    /// Set the option list.
    void setOptions(const std::vector<std::string>& opts) {
      _options = opts;
      buildOptionMap();
    }

    /// Build a map of options to facilitate checking.
    void buildOptionMap();

    /// @}


    /// @defgroup anainfo_status Status info and categories
    /// @{

    /// Whether this analysis is trusted (in any way!)
    const std::string& status() const { return _status; }
    /// Set the analysis code status.
    void setStatus(const std::string& status) { _status = status; }

    /// Return true if finalize() can be run multiple times for this analysis.
    bool reentrant() const { return _reentrant; }
    /// Set re-entrant status
    void setReentrant(bool ree=true) { _reentrant = ree; }

    /// Return true if validated
    bool validated() const {
      return statuscheck("VALIDATED");
    }

    /// Return true if preliminary
    bool preliminary() const {
      return statuscheck("PRELIMINARY");
    }

    /// Return true if obsolete
    bool obsolete() const {
      return statuscheck("OBSOLETE");
    }

    /// Return true if unvalidated
    bool unvalidated() const {
      return statuscheck("UNVALIDATED");
    }

    /// Return true if includes random variations
    bool random() const {
      return statuscheck("RANDOM");
    }

    /// Return true if the analysis uses generator-dependent information.
    bool unphysical() const {
      return statuscheck("UNPHYSICAL");
    }

    /// Check if refdata comes automatically from Hepdata.
    bool hepdata() const {
      return !statuscheck("NOHEPDATA");
    }

    /// Check if this analysis can handle multiple weights.
    bool multiweight() const {
      return !statuscheck("SINGLEWEIGHT");
    }

    /// Helper function for checking status-string contents
    bool statuscheck(const string& word) const {
      auto pos =_status.find(word);
      if ( pos == string::npos ) return false;
      if ( pos > 0 && isalnum(_status[pos - 1]) ) return false;
      if ( pos + word.length() < _status.length() &&
           isalnum(_status[pos + word.length()]) ) return false;
      return true;
    }

    /// @}


    /// Find the path to the reference-data file for this analysis
    std::string refFile() const;

    /// List a series of command lines to be used for valdation
    const std::vector<std::string> & validation() const {
      return _validation;
    }

    /// Return true if this analysis needs to know the process cross-section.
    ///
    /// @deprecated Cross-section should now always be available from the HepMC
    bool needsCrossSection() const { return _needsCrossSection; }


  private:

    // std::map<string,string> _yamldict;

    std::string _name;
    std::string _refDataName;
    std::string _spiresId, _inspireId;
    std::vector<std::string> _authors;
    std::string _summary;
    std::string _description;
    std::string _runInfo;
    std::string _experiment;
    std::string _collider;
    std::vector<std::pair<PdgId, PdgId> > _beams;
    std::vector<std::pair<double, double> > _energies;
    std::string _year;
    double _luminosityfb;
    std::vector<std::string> _references;
    std::vector<std::string> _keywords;
    std::string _bibKey;
    std::string _bibTeX;
    //std::string _bibTeXBody; ///< Was thinking of avoiding duplication of BibKey...
    std::string _status;
    std::string _warning;
    std::string _refmatch;
    std::string _refunmatch;
    std::string _writerdoubleprecision;
    std::vector<std::string> _todos;
    bool _needsCrossSection;

    std::vector<std::string> _options;
    std::map< std::string, std::set<std::string> > _optionmap;

    std::vector<std::string> _validation;

    bool _reentrant;

    void clear() {
      //_yamldict.clear();
      _name = "";
      _refDataName = "";
      _spiresId = "";
      _inspireId = "";
      _authors.clear();
      _summary = "";
      _description = "";
      _runInfo = "";
      _experiment = "";
      _collider = "";
      _beams.clear();
      _energies.clear();
      _year = "";
      _luminosityfb = -1;
      _references.clear();
      _keywords.clear();
      _bibKey = "";
      _bibTeX = "";
      //_bibTeXBody = "";
      _status = "";
      _warning = "";
      _refmatch = "";
      _refunmatch = "";
      _writerdoubleprecision = "";
      _todos.clear();
      _needsCrossSection = false;
      _options.clear();
      _optionmap.clear();
      _validation.clear();
      _reentrant = false;
    }

  };


  /// String representation
  std::string toString(const AnalysisInfo& ai);

  /// Stream an AnalysisInfo as a text description
  inline std::ostream& operator<<(std::ostream& os, const AnalysisInfo& ai) {
    os << toString(ai);
    return os;
  }


}

#endif

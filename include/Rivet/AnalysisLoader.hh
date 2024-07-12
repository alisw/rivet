// -*- C++ -*-
#ifndef RIVET_AnalysisLoader_HH
#define RIVET_AnalysisLoader_HH

#include "Rivet/Config/RivetCommon.hh"
#include <map>
#include <string>

namespace Rivet {


  // Forward declarations
  class Analysis;
  class AnalysisBuilderBase;
  class Log;


  /// Internal class which loads and registers analyses from plugin libs
  class AnalysisLoader {
  public:

    /// @name Analysis names
    /// @{

    /// Get the available analyses' canonical names
    static vector<string> analysisNames();

    /// Get all the available analyses' names, including aliases
    static vector<string> allAnalysisNames();
    /// @deprecated Use allAnalysisNames()
    static vector<string> getAllAnalysisNames() { return allAnalysisNames(); }

    /// Get the standard analyses' names (from a release-specific list file)
    static vector<string> stdAnalysisNames();

    /// Get the map of analysis alias-names to their canonical equivalents
    static map<string,string> analysisNameAliases();

    /// @}


    /// @name Analysis instantiation
    /// @{

    /// @brief Get an analysis by name
    ///
    /// Warning: a name arg which matches no known analysis will return a null
    /// pointer. Check your return values before using them!
    static unique_ptr<Analysis> getAnalysis(const string& analysisname);

    /// Get all the available analyses.
    static vector<unique_ptr<Analysis>> getAllAnalyses();

    /// @}


    /// @name Plugin library management
    /// @{

    /// @brief Return the active set of analysis plugin paths
    ///
    /// If the current set of plugin paths is empty, this will automatically
    /// call the search function.
    static vector<string> analysisPlugins();

    /// @brief Explicitly trigger the search for the available analyses plugin libraries (caches)
    ///
    /// Search the analysis paths for analysis libraries with the Rivet*.so name
    /// pattern, if the RIVET_ANALYSIS_PLUGINS environment variable has not been
    /// set. If it has, use the space-separated file paths in that variable.
    static vector<string> searchAnalysisPlugins();

    /// @brief Set a fixed list of analysis plugin libraries, bypassing the search.
    ///
    /// Setting an empty list of plugin libraries will reenable plugin searching.
    ///
    /// May be useful on MPI machines, where you want to avoid having many
    /// machines all doing the same filesystem lookup: either specify all the
    /// paths in advance, or (if you have a shared or guaranteed identical
    /// filesystem) retrieve the search function on the main rank, and pass the
    /// list for rapid, search-free setting on all the others.
    ///
    /// Calling this function will clear the lists of analysis builder functions,
    /// if any have been loaded from the previous active set of plugin libraries.
    static void setAnalysisPlugins(const vector<string> pluginpaths);

    /// Load the available analyses from the active plugin libraries (caches)
    static void loadFromAnalysisPlugins();

    /// @}


  private:

    /// Allow the analysis builders to call the private _registerBuilder function
    friend class AnalysisBuilderBase;

    /// Register a new analysis builder
    static void _registerBuilder(const AnalysisBuilderBase* ab);

    /// List of Rivet*.so plugin library paths to load from
    static vector<string> _pluginpaths;

    typedef map<string, const AnalysisBuilderBase*> AnalysisBuilderMap;
    /// Canonical analysis builder functors
    static AnalysisBuilderMap _ptrs;
    /// Alias analysis builder functors
    static AnalysisBuilderMap _aliasptrs;

  };


}

#endif

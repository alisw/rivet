// -*- C++ -*-
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/osdir.hh"
#include "Rivet/Analysis.hh"
#include <fstream>
#include <dlfcn.h>

namespace Rivet {

  using namespace std;


  // Initialise static ptr collection
  AnalysisLoader::AnalysisBuilderMap AnalysisLoader::_ptrs;
  AnalysisLoader::AnalysisBuilderMap AnalysisLoader::_aliasptrs;


  // Provide a logger function
  namespace {
    inline Log& getLog() {
      return Log::getLog("Rivet.AnalysisLoader");
    }
  }


  vector<string> AnalysisLoader::analysisNames() {
    _loadAnalysisPlugins();
    vector<string> names;
    for (const AnalysisBuilderMap::value_type& p : _ptrs) {
      names += p.second->name();
      // const string cname = p.second->name();
      // if (!contains(names, cname)) names += cname; //< avoid duplicates from alias entries in the map
    }
    return names;
  }


  vector<string> AnalysisLoader::allAnalysisNames() {
    _loadAnalysisPlugins();
    vector<string> names;
    for (const AnalysisBuilderMap::value_type& p : _ptrs) names += p.first;
    for (const AnalysisBuilderMap::value_type& p : _aliasptrs) names += p.first;
    return names;
  }


  vector<string> AnalysisLoader::stdAnalysisNames() {
    vector<string> rtn;
    const string anadatpath = findAnalysisDataFile("analyses.dat");
    if (fileexists(anadatpath)) {
      ifstream anadat(anadatpath);
      string ananame;
      while (anadat >> ananame) rtn += ananame;
    }
    return rtn;
  }


  map<string,string> AnalysisLoader::analysisNameAliases() {
    _loadAnalysisPlugins();
    map<string,string> alias_names;
    for (const AnalysisBuilderMap::value_type& p : _ptrs) {
      const string alias = p.second->alias();
      if (alias != "") alias_names[alias] = p.second->name();
    }
    return alias_names;
  }



  unique_ptr<Analysis> AnalysisLoader::getAnalysis(const string& analysisname) {
    _loadAnalysisPlugins();
    AnalysisBuilderMap::const_iterator ai = _ptrs.find(analysisname);
    if (ai == _ptrs.end()) {
      ai = _aliasptrs.find(analysisname);
      if (ai == _aliasptrs.end()) return nullptr;
      MSG_WARNING("Instantiating analysis '" << ai->second->name() << "' via alias '"
                  << analysisname << "'. Using the canonical name is recommended");
    }
    return ai->second->mkAnalysis();
  }


  vector<unique_ptr<Analysis>> AnalysisLoader::getAllAnalyses() {
    _loadAnalysisPlugins();
    vector<unique_ptr<Analysis>> analyses;
    for (const auto & p : _ptrs) {
      analyses.emplace_back( p.second->mkAnalysis() );
    }
    return analyses;
  }


  void AnalysisLoader::_registerBuilder(const AnalysisBuilderBase* ab) {
    if (!ab) return;

    // Register by canonical name
    const string name = ab->name();
    if (_ptrs.find(name) != _ptrs.end()) {
      // Duplicate analyses will be ignored... loudly
      //cerr << "Ignoring duplicate plugin analysis called '" << name << "'" << '\n';
      MSG_WARNING("Ignoring duplicate plugin analysis called '" << name << "'");
    } else {
      MSG_TRACE("Registering a plugin analysis called '" << name << "'");
      _ptrs[name] = ab;
    }

    // Register by alias name
    const string aname = ab->alias();
    if (!aname.empty()) {
      //MSG_WARNING("ALIAS!!! " << aname);
      if (_ptrs.find(aname) != _ptrs.end()) {
        MSG_WARNING("Ignoring duplicate plugin analysis alias '" << aname << "'");
      } else if (_aliasptrs.find(aname) != _aliasptrs.end()) {
        MSG_WARNING("Ignoring duplicate plugin analysis alias '" << aname << "'");
      } else {
        MSG_TRACE("Registering a plugin analysis via alias '" << aname << "'");
        _aliasptrs[aname] = ab;
      }
    }
  }


  void AnalysisLoader::_loadAnalysisPlugins() {
    // Only run once
    if (!_ptrs.empty()) return;

    // Build the list of directories to search
    const vector<string> dirs = getAnalysisLibPaths();

    // Find plugin module library files
    const string libsuffix = ".so";
    vector<string> pluginfiles;
    for (const string& d : dirs) {
      if (d.empty()) continue;
      oslink::directory dir(d);
      while (dir) {
        string filename = dir.next();
        // Require that plugin lib name starts with 'Rivet'
        if (filename.find("Rivet") != 0) continue;
        size_t posn = filename.find(libsuffix);
        if (posn == string::npos || posn != filename.length()-libsuffix.length()) continue;
        /// @todo Make sure this is an abs path
        /// @todo Sys-dependent path separator instead of "/"
        const string path = d + "/" + filename;
        // Ensure no duplicate paths
        if (find(pluginfiles.begin(), pluginfiles.end(), path) == pluginfiles.end()) {
          pluginfiles += path;
        }
      }
    }

    // Load the plugin files
    MSG_TRACE("Candidate analysis plugin libs: " << pluginfiles);
    for (const string& pf : pluginfiles) {
      MSG_TRACE("Trying to load plugin analyses from file " << pf);
      void* handle = dlopen(pf.c_str(), RTLD_LAZY);
      if (!handle) {
        MSG_WARNING("Cannot open " << pf << ": " << dlerror());
        continue;
      }
    }
  }


}

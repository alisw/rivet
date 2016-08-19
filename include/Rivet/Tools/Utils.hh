// -*- C++ -*-
#ifndef RIVET_Utils_HH
#define RIVET_Utils_HH

#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/PrettyPrint.hh"
#include <sstream>
#include <cctype>
#include <algorithm>
#include <cerrno>

namespace Rivet {


  /// @name String utils
  //@{

  struct bad_lexical_cast : public std::runtime_error {
    bad_lexical_cast(const std::string& what) : std::runtime_error(what) {}
  };

  /// @brief Convert between any types via stringstream
  template<typename T, typename U>
  T lexical_cast(const U& in) {
    try {
      std::stringstream ss;
      ss << in;
      T out;
      ss >> out;
      return out;
    } catch (const std::exception& e) {
      throw bad_lexical_cast(e.what());
    }
  }

  /// @brief Convert any object to a string
  ///
  /// Just a convenience wrapper for the more general Boost lexical_cast
  template <typename T>
  inline string to_str(const T& x) {
    return lexical_cast<string>(x);
  }

  /// @brief Convert any object to a string
  ///
  /// An alias for to_str() with a more "Rivety" mixedCase name.
  template <typename T>
  inline string toString(const T& x) {
    return to_str(x);
  }

  /// Replace the first instance of patt with repl
  inline string& replace_first(string& str, const string& patt, const string& repl) {
    if (!contains(str, patt)) return str; //< contains from RivetSTL
    str.replace(str.find(patt), patt.size(), repl);
    return str;
  }

  /// @brief Replace all instances of patt with repl
  ///
  /// @note Finding is interleaved with replacement, so the second search happens after
  /// first replacement, etc. This could lead to infinite loops and other counterintuitive
  /// behaviours if not careful.
  inline string& replace_all(string& str, const string& patt, const string& repl) {
    if (!contains(str, patt)) return str; //< contains from RivetSTL
    while (true) {
      string::size_type it = str.find(patt);
      if (it == string::npos) break;
      str.replace(it, patt.size(), repl);
    }
    return str;
  }


  /// Case-insensitive string comparison function
  inline int nocase_cmp(const string& s1, const string& s2) {
    string::const_iterator it1 = s1.begin();
    string::const_iterator it2 = s2.begin();
    while ( (it1 != s1.end()) && (it2 != s2.end()) ) {
      if(::toupper(*it1) != ::toupper(*it2)) { // < Letters differ?
        // Return -1 to indicate smaller than, 1 otherwise
        return (::toupper(*it1) < ::toupper(*it2)) ? -1 : 1;
      }
      // Proceed to the next character in each string
      ++it1;
      ++it2;
    }
    size_t size1 = s1.size(), size2 = s2.size(); // Cache lengths
    // Return -1,0 or 1 according to strings' lengths
    if (size1 == size2) return 0;
    return (size1 < size2) ? -1 : 1;
  }


  /// Case-insensitive string equality function
  inline bool nocase_equals(const string& s1, const string& s2) {
    return nocase_cmp(s1, s2) == 0;
  }


  /// Convert a string to lower-case
  inline string toLower(const string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) tolower);
    return out;
  }


  /// Convert a string to upper-case
  inline string toUpper(const string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) toupper);
    return out;
  }


  /// Check whether a string @a start is found at the start of @a s
  inline bool startsWith(const string& s, const string& start) {
    if (s.length() < start.length()) return false;
    return s.substr(0, start.length()) == start;
  }


  /// Check whether a string @a end is found at the end of @a s
  inline bool endsWith(const string& s, const string& end) {
    if (s.length() < end.length()) return false;
    return s.substr(s.length() - end.length()) == end;
  }


  /// Make a string containing the string representations of each item in v, separated by sep
  template <typename T>
  inline string join(const vector<T>& v, const string& sep=" ") {
    string rtn;
    for (size_t i = 0; i < v.size(); ++i) {
      if (i != 0) rtn += sep;
      rtn += to_str(v[i]);
    }
    return rtn;
  }

  /// Make a string containing the string representations of each item in s, separated by sep
  template <typename T>
  inline string join(const set<T>& s, const string& sep=" ") {
    string rtn;
    for (const T& x : s) {
      if (rtn.size() > 0) rtn += sep;
      rtn += to_str(x);
    }
    return rtn;
  }

  //@}


  /// @name Path utils
  //@{

  /// @brief Split a path string with colon delimiters
  ///
  /// Ignores zero-length substrings. Designed for getting elements of filesystem paths, naturally.
  inline vector<string> pathsplit(const string& path) {
    const string delim = ":";
    vector<string> dirs;
    string tmppath = path;
    while (true) {
      const size_t delim_pos = tmppath.find(delim);
      if (delim_pos == string::npos) break;
      const string dir = tmppath.substr(0, delim_pos);
      if (dir.length()) dirs.push_back(dir); // Don't insert "empties"
      tmppath.replace(0, delim_pos+1, "");
    }
    if (tmppath.length()) dirs.push_back(tmppath); // Don't forget the trailing component!
    return dirs;
  }


  /// @brief Join several filesystem paths together with the standard ':' delimiter
  ///
  /// Note that this does NOT join path elements together with a platform-portable
  /// directory delimiter, cf. the Python @c {os.path.join}!
  inline string pathjoin(const vector<string>& paths) {
    return join(paths, ":");
  }

  //@}


  /// @name Container utils
  //@{

  /// Return true if f(x) is true for any x in container c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool any(const CONTAINER& c, const FN& f) {
    for (const auto& x : c)
      if (f(x)) return true;
    return false;
  }

  /// Return true if @a f(x) is true for all @c x in container @a c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool all(const CONTAINER& c, const FN& f) {
    for (const auto& x : c)
      if (!f(x)) return false;
    return true;
  }

  /// Generic sum function, adding @a fn(@c x) for all @c x in container @a c, starting with @a start
  template <typename CONTAINER, typename FN, typename T>
  inline T sum(const CONTAINER& c, const FN& f, const T& start=T()) {
    T rtn = start;
    for (const auto& x : c)
      rtn += f(x);
    return rtn;
  }

  //@}


}

#endif

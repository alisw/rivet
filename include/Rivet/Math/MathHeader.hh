#ifndef RIVET_Math_MathHeader
#define RIVET_Math_MathHeader

#include "Rivet/Exceptions.hh"
#include <stdexcept>
#include <string>
#include <ostream>
#include <sstream>
#include <iostream>
#include <limits>
#include <climits>
#include <cfloat>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>


// Macro to help with overzealous compiler warnings
/// @note It's easier and better to just not give an arg name to args which won't be used, when possible.
#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#else
# define UNUSED(x) x
#endif


/// Macro to help mark code as deprecated to produce compiler warnings
#ifndef DEPRECATED
#if __GNUC__ && __cplusplus && RIVET_NO_DEPRECATION_WARNINGS == 0
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION >= 40500
  #if __cplusplus > 201103L
  #define DEPRECATED(x) [[deprecated(x)]]
  #else
  #define DEPRECATED(x) __attribute__((deprecated(x)))
  #endif
#else
  #define DEPRECATED(x) __attribute__((deprecated))
#endif
#else
  #define DEPRECATED(x)
#endif
#endif


namespace Rivet {

  using std::string;
  using std::ostream;
  using std::ostringstream;
  using std::cout;
  using std::endl;
  using std::pair;
  using std::vector;
  using std::transform;
  using std::min;
  using std::max;
  using std::abs;
  using std::isnan;
  using std::isinf;

  /// Pre-defined numeric type limits
  /// @deprecated Prefer the standard DBL/INT_MAX
  static const double MAXDOUBLE = DBL_MAX; // was std::numeric_limits<double>::max(); -- warns in GCC5
  static const double MAXINT = INT_MAX; // was std::numeric_limits<int>::max(); -- warns in GCC5

  /// A pre-defined value of \f$ \pi \f$.
  static const double PI = M_PI;

  /// A pre-defined value of \f$ 2\pi \f$.
  static const double TWOPI = 2*M_PI;

  /// A pre-defined value of \f$ \pi/2 \f$.
  static const double HALFPI = M_PI_2;

  /// Enum for signs of numbers.
  enum Sign { MINUS = -1, ZERO = 0, PLUS = 1 };

  /// Enum for rapidity variable to be used in calculating \f$ R \f$, applying rapidity cuts, etc.
  enum RapScheme { PSEUDORAPIDITY = 0, ETARAP = 0, RAPIDITY = 1, YRAP = 1 };

  /// Enum for range of \f$ \phi \f$ to be mapped into
  enum PhiMapping { MINUSPI_PLUSPI, ZERO_2PI, ZERO_PI };

}

#endif

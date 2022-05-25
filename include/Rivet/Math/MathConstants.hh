#ifndef RIVET_Math_MathConstants
#define RIVET_Math_MathConstants

#include "Rivet/Tools/Exceptions.hh"
#include "Rivet/Tools/Utils.hh"
#include <cmath>

namespace Rivet {


  /// Pre-defined numeric type limits
  /// A pre-defined value of \f$ \pi \f$.
  static const double PI = M_PI;

  /// A pre-defined value of \f$ 2\pi \f$.
  static const double TWOPI = 2*M_PI;

  /// A pre-defined value of \f$ \pi/2 \f$.
  static const double HALFPI = M_PI_2;

  /// A pre-defined value of \f$ \sqrt{2} \f$.
  static const double SQRT2 = M_SQRT2;

  /// A pre-defined value of \f$ \sqrt{\pi} \f$.
  static const double SQRTPI = 2 / M_2_SQRTPI;

  // /// A pre-defined value of \f$ \sqrt{2\pi} \f$.
  // static const double SQRT2PI = SQRT2 * SQRTPI;

  /// @brief Pre-defined values of \f$ \infty \f$.
  ///
  /// See https://en.cppreference.com/w/cpp/types/numeric_limits/infinity
  static const double INFF = HUGE_VALF;
  static const double INF = HUGE_VAL;
  static const double INFL = HUGE_VALL;

  // Other useful predefined values already exist in C++, e.g.:
  // DBL_MAX
  // NAN


  /// Enum for signs of numbers.
  enum Sign { MINUS = -1, ZERO = 0, PLUS = 1 };

  /// Enum for rapidity variable to be used in calculating \f$ R \f$, applying rapidity cuts, etc.
  enum RapScheme { PSEUDORAPIDITY = 0, ETARAP = 0, RAPIDITY = 1, YRAP = 1 };

  /// Enum for range of \f$ \phi \f$ to be mapped into
  enum PhiMapping { MINUSPI_PLUSPI, ZERO_2PI, ZERO_PI };

}

#endif

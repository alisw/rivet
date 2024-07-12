#ifndef RIVET_MATH_UNITS
#define RIVET_MATH_UNITS

#include "Rivet/Math/MathConstants.hh"

namespace Rivet {

  //
  // Length [L]
  //
  constexpr double millimeter  = 1.;
  constexpr double millimeter2 = millimeter*millimeter;
  constexpr double millimeter3 = millimeter*millimeter*millimeter;

  constexpr double centimeter  = 10.*millimeter;
  constexpr double centimeter2 = centimeter*centimeter;
  constexpr double centimeter3 = centimeter*centimeter*centimeter;

  constexpr double meter  = 1000.*millimeter;
  constexpr double meter2 = meter*meter;
  constexpr double meter3 = meter*meter*meter;

  constexpr double micrometer = 1.e-6 *meter;
  constexpr double nanometer  = 1.e-9 *meter;
  constexpr double angstrom   = 1.e-10*meter;
  constexpr double picometer  = 1.e-12*meter;
  constexpr double femtometer = 1.e-15*meter;
  constexpr double attometer  = 1.e-18*meter;
  constexpr double fermi      = femtometer;

  // symbols
  constexpr double mm  = millimeter;
  constexpr double mm2 = millimeter2;
  constexpr double mm3 = millimeter3;

  constexpr double cm  = centimeter;
  constexpr double cm2 = centimeter2;
  constexpr double cm3 = centimeter3;

  constexpr double m  = meter;
  constexpr double m2 = meter2;
  constexpr double m3 = meter3;

  // constexpr double barn = 1.e-28*meter2;
  // Barn-units in terms of the pb returned by AGILe
  constexpr double  picobarn = 1.0;
  constexpr double      barn = 1.0e+12* picobarn;
  constexpr double millibarn = 1.0e-3 * barn;
  constexpr double microbarn = 1.0e-6 * barn;
  constexpr double  nanobarn = 1.0e-9 * barn;
  constexpr double femtobarn = 1.0e-15 * barn;
  constexpr double attobarn  = 1.0e-18 * barn;

  //
  // Time [T]
  //
  constexpr double nanosecond  = 1.0;
  constexpr double second      = 1.e+9 *nanosecond;
  constexpr double millisecond = 1.e-3 *second;
  constexpr double microsecond = 1.e-6 *second;
  constexpr double  picosecond = 1.e-12*second;

  // symbols
  constexpr double ns = nanosecond;
  constexpr double  s = second;
  constexpr double ms = millisecond;

  //
  // Electric charge [Q]
  //
  constexpr double eplus = 1.0;		// positron charge
  constexpr double e_SI  = 1.60217733e-19;	// positron charge in coulomb

  //
  // Energy [E]
  //
  constexpr double gigaelectronvolt = 1.;
  constexpr double     electronvolt = 1.e-9*gigaelectronvolt;
  constexpr double kiloelectronvolt = 1.e-6*gigaelectronvolt;
  constexpr double megaelectronvolt = 1.e-3*gigaelectronvolt;
  constexpr double teraelectronvolt = 1.e+3*gigaelectronvolt;
  constexpr double petaelectronvolt = 1.e+6*gigaelectronvolt;

  // symbols
  constexpr double  eV = electronvolt;
  constexpr double keV = kiloelectronvolt;
  constexpr double MeV = megaelectronvolt;
  constexpr double GeV = gigaelectronvolt;
  constexpr double TeV = teraelectronvolt;
  constexpr double PeV = petaelectronvolt;

  constexpr double  eV2 = eV*eV;
  constexpr double keV2 = keV*keV;
  constexpr double MeV2 = MeV*MeV;
  constexpr double GeV2 = GeV*GeV;
  constexpr double TeV2 = TeV*TeV;
  constexpr double PeV2 = PeV*PeV;

}

#endif

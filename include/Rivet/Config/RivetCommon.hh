#ifndef RIVET_RivetCommon_HH
#define RIVET_RivetCommon_HH

// Convenience build-setup header for Rivet internal use


// Automatic build info from autoconf
#include "Rivet/Config/RivetConfig.hh"
#include "Rivet/Config/BuildOptions.hh"


/// Macro to help with overzealous compiler warnings
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


#include "Rivet/Exceptions.hh"

#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/RivetBoost.hh"
#include "Rivet/Tools/RivetHepMC.hh"

#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/Logging.hh"

#include "Rivet/Math/Units.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Constants.hh"

// #include "Rivet/Particle.hh"
// #include "Rivet/Event.hh"

#endif

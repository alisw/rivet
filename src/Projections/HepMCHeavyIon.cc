// -*- C++ -*-
#include "Rivet/Projections/HepMCHeavyIon.hh"

#ifdef RIVET_ENABLE_HEPMC_3
#define IMPLEMENTATION(rettype, functionname, defret) \
rettype HepMCHeavyIon::functionname() const { \
  return _hi? _hi->functionname: defret;        \
}

#define IMPLEMENTATION_ALT_HEPMC2(rettype, functionname, hepmc2name, defret) \
rettype HepMCHeavyIon::functionname() const { \
  return _hi? _hi->functionname: defret;        \
}

#define IMPLEMENTATION_NO_HEPMC2(rettype, functionname, defret) \
rettype HepMCHeavyIon::functionname() const { \
  return _hi? _hi->functionname: defret;        \
}

#else

#define IMPLEMENTATION(rettype, functionname, defret) \
rettype HepMCHeavyIon::functionname() const { \
  return _hi? _hi->functionname(): defret;      \
}

#define IMPLEMENTATION_ALT_HEPMC2(rettype, functionname, hepmc2name, defret) \
rettype HepMCHeavyIon::functionname() const { \
  return _hi? _hi->, hepmc2name(): defret;      \
}

#define IMPLEMENTATION_NO_HEPMC2(rettype, functionname, defret) \
rettype HepMCHeavyIon::functionname() const { \
 MSG_WARNING("HeavyIon::" #functionname " is only avialable in HepMC3"); \
 return defret; \
}

#endif


namespace Rivet {


HepMCHeavyIon::HepMCHeavyIon() {
  setName("HepMCHeavyIon");
}

void HepMCHeavyIon::project(const Event& e) {
  _hi = e.genEvent()->heavy_ion();
  if ( !_hi )
    MSG_WARNING("Could not find the HepMC HeavyIon object");
}

IMPLEMENTATION(int, Ncoll_hard, -1)

IMPLEMENTATION(int, Npart_proj, -1)

IMPLEMENTATION(int, Npart_targ, -1)

IMPLEMENTATION(int, Ncoll, -1)

IMPLEMENTATION(int, N_Nwounded_collisions, -1)

IMPLEMENTATION(int, Nwounded_N_collisions, -1)

IMPLEMENTATION(int, Nwounded_Nwounded_collisions, -1)

IMPLEMENTATION(double, impact_parameter, -1.0)

IMPLEMENTATION(double, event_plane_angle, -1.0)

IMPLEMENTATION(double, sigma_inel_NN, -1.0)

#ifdef RIVET_ENABLE_HEPMC_20610
IMPLEMENTATION(double, centrality, -1.0)
#else
IMPLEMENTATION_NO_HEPMC2(double, centrality, -1.0)
#endif

IMPLEMENTATION_NO_HEPMC2(double, user_cent_estimate, -1.0)

IMPLEMENTATION_NO_HEPMC2(int, Nspec_proj_n, -1)

IMPLEMENTATION_NO_HEPMC2(int, Nspec_targ_n, -1)

IMPLEMENTATION_NO_HEPMC2(int, Nspec_proj_p, -1)

IMPLEMENTATION_NO_HEPMC2(int, Nspec_targ_p, -1)

map<int,double> HepMCHeavyIon::participant_plane_angles() const {
#ifdef RIVET_ENABLE_HEPMC_3
  return _hi? _hi->participant_plane_angles: map<int,double>(); 
#else
  MSG_WARNING("HeavyIon::participant_plane_angles is only avialable in HepMC3");
  return map<int,double>(); 
#endif
}

map<int,double> HepMCHeavyIon::eccentricities() const {
#ifdef RIVET_ENABLE_HEPMC_3
  return _hi? _hi->eccentricities: map<int,double>(); 
#else
  MSG_WARNING("HeavyIon::eccentricities is only avialable in HepMC3");
  return map<int,double>(); 
#endif
}

}

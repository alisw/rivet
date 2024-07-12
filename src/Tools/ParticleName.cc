#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/Utils.hh"

namespace Rivet {

  namespace PID {


    // Initialise ParticleNames singleton pointer
    unique_ptr<ParticleNames> ParticleNames::_instance = nullptr;


    std::string ParticleNames::_particleName(PdgId pid) {
      if (_ids_names.find(pid) == _ids_names.end()) {
        //throw PidError("Particle ID '" + lexical_cast<string>(pid) + "' not known.");
        return lexical_cast<string>(pid);
      }
      return _ids_names[pid];
    }


    PdgId ParticleNames::_particleId(const std::string& pname) {
      if (_names_ids.find(pname) == _names_ids.end()) {
        if (toUpper(pname) == "P+" || toUpper(pname) == "P") return PROTON;
        if (toUpper(pname) == "P-" || toUpper(pname) == "PBAR") return ANTIPROTON;
        if (toUpper(pname) == "E-") return ELECTRON;
        if (toUpper(pname) == "E+") return POSITRON;
        if (toUpper(pname) == "GAMMA") return PHOTON;
        if (toUpper(pname) == "N") return NEUTRON;
        if (toUpper(pname) == "D") return DEUTERON;
        if (toUpper(pname) == "AL") return ALUMINIUM;
        if (toUpper(pname) == "CU") return COPPER;
        if (toUpper(pname) == "XE") return XENON;
        if (toUpper(pname) == "AU") return GOLD;
        if (toUpper(pname) == "PB") return LEAD;
        if (toUpper(pname) == "U") return URANIUM;
        try {
          PdgId rtn = lexical_cast<PdgId>(pname);
          return rtn;
        } catch (const bad_lexical_cast& blc) {
          //throw PidError("Particle name '" + pname + "' not known and could not be directly cast to a PDG ID.");
          return 0;
        }
      }
      return _names_ids[pname];
    }


  }

}

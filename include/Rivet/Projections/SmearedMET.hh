// -*- C++ -*-
#ifndef RIVET_SmearedMET_HH
#define RIVET_SmearedMET_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  /// Wrapper projection for smearing missing (transverse) energy/momentum with detector resolutions
  class SmearedMET : public Projection {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor from a MissingMomentum projection and a smearing function
    template <typename V2VFN>
    SmearedMET(const MissingMomentum& mm, const V2VFN& metSmearFn)
      : _metSmearFn(metSmearFn)
    {
      setName("SmearedMET");
      addProjection(mm, "TruthMET");
    }

    /// @brief Constructor from a Cut (on the particles used to determine missing momentum) and a smearing function
    template <typename V2VFN>
    SmearedMET(const V2VFN& metSmearFn, const Cut& cut)
      : _metSmearFn(metSmearFn)
    {
      setName("SmearedMET");
      addProjection(MissingMomentum(cut), "TruthMET");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedMET);

    //@}


    /// Compare to another SmearedMET
    int compare(const Projection& p) const {
      const SmearedMET& other = dynamic_cast<const SmearedMET&>(p);
      if (get_address(_metSmearFn) == 0) return UNDEFINED;
      MSG_TRACE("Smear hashes = " << get_address(_metSmearFn) << "," << get_address(other._metSmearFn));
      return mkPCmp(other, "TruthMET") || cmp(get_address(_metSmearFn), get_address(other._metSmearFn));
    }


    /// Perform the MET finding & smearing calculation
    void project(const Event& e) {
      _vet = applyProjection<MissingMomentum>(e, "TruthMET").vectorEt();
      if (_metSmearFn) _vet = _metSmearFn(_vet); //< smearing
    }


    /// The vector-summed visible transverse energy in the event, as a 3-vector with z=0
    /// @note Reverse this vector with operator- to get the missing ET vector.
    const Vector3& vectorEt() const { return _vet; }

    /// The vector-summed missing transverse energy in the event.
    double missingEt() const { return vectorEt().mod(); }
    /// Alias for missingEt
    double met() const { return missingEt(); }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() {  }


  private:

    Vector3 _vet;

    /// Stored smearing function
    std::function<Vector3(const Vector3&)> _metSmearFn;

  };


}

#endif

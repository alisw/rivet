// -*- C++ -*-
#ifndef RIVET_SmearedMET_HH
#define RIVET_SmearedMET_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/METFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  /// Wrapper projection for smearing missing (transverse) energy/momentum with detector resolutions
  class SmearedMET : public METFinder {
  public:

    /// @name Constructors etc.
    ///@{

    /// @brief Constructor from a MissingMomentum projection and a smearing function
    ///
    /// Smearing function maps a 3-vector MET and scalar SET to a new MET 3-vector: f(V3, double) -> V3
    template <typename V2VFN>
    SmearedMET(const MissingMomentum& mm, const V2VFN& metSmearFn)
      : _metSmearFn(metSmearFn)
    {
      setName("SmearedMET");
      declare(mm, "TruthMET");
    }

    /// @brief Constructor from a Cut (on the particles used to determine missing momentum) and a smearing function
    template <typename V2VFN>
    SmearedMET(const V2VFN& metSmearFn, const Cut& cut)
      : _metSmearFn(metSmearFn)
    {
      setName("SmearedMET");
      declare(MissingMomentum(cut), "TruthMET");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedMET);

    ///@}


    /// Compare to another SmearedMET
    CmpState compare(const Projection& p) const {
      const SmearedMET& other = dynamic_cast<const SmearedMET&>(p);
      if (get_address(_metSmearFn) == 0) return cmp((size_t)this, (size_t)&p);
      MSG_TRACE("Smear hashes = " << get_address(_metSmearFn) << "," << get_address(other._metSmearFn));
      return mkPCmp(other, "TruthMET") || cmp(get_address(_metSmearFn), get_address(other._metSmearFn));
    }


    /// Perform the MET finding & smearing calculation
    void project(const Event& e) {
      const auto& mm = apply<MissingMomentum>(e, "TruthMET");
      _vet = mm.vectorEt();
      if (_metSmearFn) _vet = _metSmearFn(_vet, mm.scalarEt()); //< smearing
    }


    /// @name Transverse momentum functions
    ///
    /// @note This may be what you want, even if the paper calls it "missing Et"!
    ///@{

    /// The vector-summed visible transverse momentum in the event, as a 3-vector with z=0
    ///
    /// @note Reverse this vector with operator- to get the missing pT vector.
    /// @todo Currently equivalent to vectorEt
    const Vector3& vectorPt() const { return vectorEt(); }

    ///@}


    /// @name Transverse energy functions
    ///
    /// @warning Despite the common names "MET" and "SET", what's often meant is the pT functions above!
    ///@{

    /// The vector-summed visible transverse energy in the event, as a 3-vector with z=0
    ///
    /// @note Reverse this vector with operator- to get the missing ET vector.
    const Vector3& vectorEt() const { return _vet; }

    //@}


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() {  }


  private:

    Vector3 _vet;

    /// Stored smearing function
    std::function<Vector3(const Vector3&, double)> _metSmearFn;

  };


}

#endif

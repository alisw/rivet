// -*- C++ -*-
#ifndef RIVET_SmearedJets_HH
#define RIVET_SmearedJets_HH

#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  /// @todo Chain constructors?
  class SmearedJets : public JetAlg {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn)
      : SmearedJets(ja, jetSmearFn, JET_BTAG_PERFECT, JET_CTAG_PERFECT, JET_EFF_ONE)
    {    }


    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN, typename J2DFN>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn,
                const J2DFN& bTagEffFn)
      : SmearedJets(ja, jetSmearFn, bTagEffFn, JET_CTAG_PERFECT, JET_EFF_ONE)
    {    }


    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN, typename J2DFNa, typename J2DFNb>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn,
                const J2DFNa& bTagEffFn,
                const J2DFNb& cTagEffFn)
      : SmearedJets(ja, jetSmearFn, bTagEffFn, cTagEffFn, JET_EFF_ONE)
    {    }


    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2JFN, typename J2DFNa, typename J2DFNb, typename J2DFNc>
    SmearedJets(const JetAlg& ja,
                const J2JFN& jetSmearFn,
                const J2DFNa& bTagEffFn,
                const J2DFNb& cTagEffFn,
                const J2DFNc& jetEffFn)
      : _jetEffFn(jetEffFn), _bTagEffFn(bTagEffFn), _cTagEffFn(cTagEffFn), _jetSmearFn(jetSmearFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedJets);

    //@}


    /// Compare to another SmearedJets
    int compare(const Projection& p) const {
      const SmearedJets& other = dynamic_cast<const SmearedJets&>(p);
      if (get_address(_jetEffFn) == 0) return UNDEFINED;
      if (get_address(_bTagEffFn) == 0) return UNDEFINED;
      if (get_address(_cTagEffFn) == 0) return UNDEFINED;
      if (get_address(_jetSmearFn) == 0) return UNDEFINED;
      MSG_TRACE("Eff hashes = " << get_address(_jetEffFn) << "," << get_address(other._jetEffFn) << "; " <<
                "smear hashes = " << get_address(_jetSmearFn) << "," << get_address(other._jetSmearFn) << "; " <<
                "b-tag hashes = " << get_address(_bTagEffFn) << "," << get_address(other._bTagEffFn) << "; " <<
                "c-tag hashes = " << get_address(_cTagEffFn) << "," << get_address(other._cTagEffFn));
      return
        mkPCmp(other, "TruthJets") ||
        cmp(get_address(_jetEffFn), get_address(other._jetEffFn)) ||
        cmp(get_address(_jetSmearFn), get_address(other._jetSmearFn)) ||
        cmp(get_address(_bTagEffFn), get_address(other._bTagEffFn)) ||
        cmp(get_address(_cTagEffFn), get_address(other._cTagEffFn));
    }


    /// Perform the jet finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Jets& truthjets = applyProjection<JetAlg>(e, "TruthJets").jetsByPt();
      _recojets.clear(); _recojets.reserve(truthjets.size());
      for (const Jet& j : truthjets) {
        // Efficiency sampling
        const double jeff = _jetEffFn ? _jetEffFn(j) : 1;
        MSG_DEBUG("Efficiency of jet " << j.mom() << " = " << 100*jeff << "%");
        MSG_DEBUG("Efficiency of jet with mom=" << j.mom()/GeV << " GeV, "
                  << "pT=" << j.pT()/GeV << ", eta=" << j.eta() << " : " << 100*jeff << "%");
        if (jeff <= 0) continue; //< no need to roll expensive dice (and we deal with -ve probabilities, just in case)
        if (jeff < 1 && rand01() > jeff) continue; //< roll dice (and deal with >1 probabilities, just in case)
        // Kinematic smearing
        Jet sj = _jetSmearFn ? _jetSmearFn(j) : j;
        MSG_DEBUG("Jet smearing from " << j.mom() << " to " << sj.mom());
        // Re-add constituents & tags if (we assume accidentally) they were lost by the smearing function
        if (sj.particles().empty() && !j.particles().empty()) sj.particles() = j.particles();
        if (sj.tags().empty() && !j.tags().empty()) sj.tags() = j.tags();
        _recojets.push_back(sj);
      }
      // Apply tagging efficiencies, using smeared kinematics as input to the tag eff functions
      for (Jet& j : _recojets) {
        const double beff = _bTagEffFn ? _bTagEffFn(j) : 1;
        const bool btag = beff == 1 || (beff != 0 && beff < rand01());
        // Remove b-tags if needed, and add a dummy one if needed
        if (!btag && j.bTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasBottom), j.tags().end());
        if (btag && !j.bTagged()) j.tags().push_back(Particle(PID::BQUARK, j.mom()));
        const double ceff = _cTagEffFn ? _cTagEffFn(j) : 1;
        const bool ctag = ceff == 1 || (ceff != 0 && ceff < rand01());
        // Remove c-tags if needed, and add a dummy one if needed
        if (!ctag && j.cTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasCharm), j.tags().end());
        if (ctag && !j.cTagged()) j.tags().push_back(Particle(PID::CQUARK, j.mom()));
      }
    }


    /// Return the full jet list for the JetAlg methods to use
    Jets _jets() const { return _recojets; }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _recojets.clear(); }


  private:

    Jets _recojets;

    /// Stored efficiency functions
    std::function<double(const Jet&)> _jetEffFn, _bTagEffFn, _cTagEffFn;

    /// Stored smearing function
    std::function<Jet(const Jet&)> _jetSmearFn;

  };


}

#endif

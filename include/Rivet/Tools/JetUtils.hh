#ifndef RIVET_JETUTILS_HH
#define RIVET_JETUTILS_HH

#include "Rivet/Jet.hh"
#include "Rivet/Tools/ParticleBaseUtils.hh"

namespace Rivet {


  /// @defgroup jetutils Functions for Jets
  /// @{

  /// @defgroup jetutils_conv Converting between Jets, Particles and PseudoJets
  /// @{

  inline PseudoJets mkPseudoJets(const Particles& ps) {
    PseudoJets rtn; rtn.reserve(ps.size());
    for (const Particle& p : ps)
      rtn.push_back(p);
    return rtn;
  }

  inline PseudoJets mkPseudoJets(const Jets& js) {
    PseudoJets rtn; rtn.reserve(js.size());
    for (const Jet& j : js)
      rtn.push_back(j);
    return rtn;
  }

  inline Jets mkJets(const PseudoJets& pjs) {
    Jets rtn; rtn.reserve(pjs.size());
    for (const PseudoJet& pj : pjs)
      rtn.push_back(pj);
    return rtn;
  }

  /// @}


  /// @defgroup jetutils_j2bool Jet classifier -> bool functors
  /// @{

  /// std::function instantiation for functors taking a Jet and returning a bool
  using JetSelector = function<bool(const Jet&)>;
  /// std::function instantiation for functors taking two Jets and returning a bool
  using JetSorter = function<bool(const Jet&, const Jet&)>;


  /// Base type for Jet -> bool functors
  struct BoolJetFunctor {
    virtual bool operator()(const Jet& p) const = 0;
    virtual ~BoolJetFunctor() {}
  };


  /// Functor for and-combination of selector logic
  struct BoolJetAND : public BoolJetFunctor {
    BoolJetAND(const std::vector<JetSelector>& sels) : selectors(sels) {}
    BoolJetAND(const JetSelector& a, const JetSelector& b) : selectors({a,b}) {}
    BoolJetAND(const JetSelector& a, const JetSelector& b, const JetSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Jet& j) const {
      for (const JetSelector& sel : selectors) if (!sel(j)) return false;
      return true;
    }
    std::vector<JetSelector> selectors;
  };
  /// Operator syntactic sugar for AND construction
  inline BoolJetAND operator && (const JetSelector& a, const JetSelector& b) {
    return BoolJetAND(a, b);
  }


  /// Functor for or-combination of selector logic
  struct BoolJetOR : public BoolJetFunctor {
    BoolJetOR(const std::vector<JetSelector>& sels) : selectors(sels) {}
    BoolJetOR(const JetSelector& a, const JetSelector& b) : selectors({a,b}) {}
    BoolJetOR(const JetSelector& a, const JetSelector& b, const JetSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Jet& j) const {
      for (const JetSelector& sel : selectors) if (sel(j)) return true;
      return false;
    }
    std::vector<JetSelector> selectors;
  };
  /// Operator syntactic sugar for OR construction
  inline BoolJetOR operator || (const JetSelector& a, const JetSelector& b) {
    return BoolJetOR(a, b);
  }


  /// Functor for inverting selector logic
  struct BoolJetNOT : public BoolJetFunctor {
    BoolJetNOT(const JetSelector& sel) : selector(sel) {}
    bool operator()(const Jet& j) const { return !selector(j); }
    JetSelector selector;
  };
  /// Operator syntactic sugar for NOT construction
  inline BoolJetNOT operator ! (const JetSelector& a) {
    return BoolJetNOT(a);
  }



  /// B-tagging functor, with a tag selection cut as the stored state
  struct HasBTag : BoolJetFunctor {
    HasBTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasBTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.bTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasBTag = HasBTag;

  /// C-tagging functor, with a tag selection cut as the stored state
  struct HasCTag : BoolJetFunctor {
    HasCTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasCTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.cTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasCTag = HasCTag;

  /// Tau-tagging functor, with a tag selection cut as the stored state
  struct HasTauTag : BoolJetFunctor {
    HasTauTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasTauTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.tauTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasTauTag = HasTauTag;

  /// Anti-B/C-tagging functor, with a tag selection cut as the stored state
  struct HasNoTag : BoolJetFunctor {
    HasNoTag(const Cut& c=Cuts::open(), bool quarktagsonly=false) : cut(c), qtagsonly(quarktagsonly) {}
    // HasNoTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return !j.bTagged(cut) && !j.cTagged(cut) && (qtagsonly || !j.tauTagged(cut)); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
    bool qtagsonly;
  };
  using hasNoTag = HasNoTag;

  /// @}


  /// @defgroup jetutils_filt Unbound functions for filtering jets
  /// @{

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  Jets& ifilter_select(Jets& jets, const Cut& c);
  /// Alias for ifilter_select
  /// @deprecated Use ifilter_select
  inline Jets& ifilterBy(Jets& jets, const Cut& c) { return ifilter_select(jets, c); }
  /// New alias for ifilter_select
  inline Jets& iselect(Jets& jets, const Cut& c) { return ifilter_select(jets, c); }


  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets filter_select(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return ifilter_select(rtn, c);
  }
  /// Alias for filter_select
  /// @deprecated Use filter_select
  inline Jets filterBy(const Jets& jets, const Cut& c) { return filter_select(jets, c); }
  /// New alias for filter_select
  inline Jets select(const Jets& jets, const Cut& c) { return filter_select(jets, c); }


  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets filter_select(const Jets& jets, const Cut& c, Jets& out) {
    out = filter_select(jets, c);
    return out;
  }
  /// Alias for filter_select
  /// @deprecated Use filter_select
  inline Jets filterBy(const Jets& jets, const Cut& c, Jets& out) { return filter_select(jets, c, out); }
  /// New alias for filter_select
  inline Jets select(const Jets& jets, const Cut& c, Jets& out) { return filter_select(jets, c, out); }



  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  Jets& ifilter_discard(Jets& jets, const Cut& c);
  /// New alias for ifilter_discard
  inline Jets& idiscard(Jets& jets, const Cut& c) { return ifilter_discard(jets, c); }


  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets filter_discard(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return ifilter_discard(rtn, c);
  }
  /// New alias for filter_discard
  inline Jets discard(const Jets& jets, const Cut& c) { return filter_discard(jets, c); }


  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets filter_discard(const Jets& jets, const Cut& c, Jets& out) {
    out = filter_discard(jets, c);
    return out;
  }
  /// New alias for filter_discard
  inline Jets discard(const Jets& jets, const Cut& c, Jets& out) { return filter_discard(jets, c, out); }

  /// @}



  /// @defgroup jetutils_coll Operations on collections of Jet
  /// @note This can't be done on generic collections of ParticleBase -- thanks, C++ :-/
  /// @{
  namespace Kin {

    inline double sumPt(const Jets& js) {
      return sum(js, pT, 0.0);
    }

    inline FourMomentum sumP4(const Jets& js) {
      return sum(js, p4, FourMomentum());
    }

    inline Vector3 sumP3(const Jets& js) {
      return sum(js, p3, Vector3());
    }

    /// @todo Min dPhi, min dR?
    /// @todo Isolation routines?

  }
  /// @}


  // Import Kin namespace into Rivet
  using namespace Kin;


  /// @}

}

#endif

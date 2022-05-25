// -*- C++ -*-
#ifndef RIVET_FastJets_HH
#define RIVET_FastJets_HH

#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/RivetFastJet.hh"

#include "fastjet/SISConePlugin.hh"
#include "fastjet/ATLASConePlugin.hh"
#include "fastjet/CMSIterativeConePlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/D0RunIIConePlugin.hh"
#include "fastjet/TrackJetPlugin.hh"
#include "fastjet/JadePlugin.hh"

#include "Rivet/Projections/PxConePlugin.hh"
#include "Rivet/Tools/TypeTraits.hh"

namespace Rivet {


  /// Project out jets found using the FastJet package jet algorithms.
  class FastJets : public JetAlg {
  public:

    /// Wrapper enum for selected FastJet jet algorithms.
    /// @todo Move to JetAlg and alias here?
    enum Algo { KT=0,
                AKT=1, ANTIKT=1,
                CA=2, CAM=2, CAMBRIDGE=2,
                SISCONE, PXCONE,
                ATLASCONE, CMSCONE,
                CDFJETCLU, CDFMIDPOINT, D0ILCONE,
                JADE, DURHAM, TRACKJET, GENKTEE ,
                KTET, ANTIKTET };


    /// @name Constructors etc.
    /// @{

    /// Constructor from a FastJet JetDefinition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             const fastjet::JetDefinition& jdef,
             JetAlg::Muons usemuons=JetAlg::Muons::ALL,
             JetAlg::Invisibles useinvis=JetAlg::Invisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : JetAlg(fsp, usemuons, useinvis), _jdef(jdef), _adef(adef)
    {
      _initBase();
    }

    /// JetDefinition-based constructor with reordered args for easier specification of jet area definition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             const fastjet::JetDefinition& jdef,
             fastjet::AreaDefinition* adef,
             JetAlg::Muons usemuons=JetAlg::Muons::ALL,
             JetAlg::Invisibles useinvis=JetAlg::Invisibles::NONE)
      : FastJets(fsp, jdef, usemuons, useinvis, adef)
    {    }

    /// Native argument constructor, using FastJet alg/scheme enums.
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter,
             JetAlg::Muons usemuons=JetAlg::Muons::ALL,
             JetAlg::Invisibles useinvis=JetAlg::Invisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : FastJets(fsp, fastjet::JetDefinition(type, rparameter, recom), usemuons, useinvis, adef)
    {    }

    /// Native argument constructor with reordered args for easier specification of jet area definition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter,
             fastjet::AreaDefinition* adef,
             JetAlg::Muons usemuons=JetAlg::Muons::ALL,
             JetAlg::Invisibles useinvis=JetAlg::Invisibles::NONE)
      : FastJets(fsp, type, recom, rparameter, usemuons, useinvis, adef)
    {    }

    /// @brief Explicitly pass in an externally-constructed plugin
    ///
    /// @warning Provided plugin and area definition pointers must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             fastjet::JetDefinition::Plugin* plugin,
             JetAlg::Muons usemuons=JetAlg::Muons::ALL,
             JetAlg::Invisibles useinvis=JetAlg::Invisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : FastJets(fsp, fastjet::JetDefinition(plugin), usemuons, useinvis, adef)
    {
      _plugin.reset(plugin);
    }

    /// @brief Explicitly pass in an externally-constructed plugin, with reordered args for easier specification of jet area definition
    ///
    /// @warning Provided plugin and area definition pointers must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             fastjet::JetDefinition::Plugin* plugin,
             fastjet::AreaDefinition* adef,
             JetAlg::Muons usemuons=JetAlg::Muons::ALL,
             JetAlg::Invisibles useinvis=JetAlg::Invisibles::NONE)
      : FastJets(fsp, plugin, usemuons, useinvis, adef)
    {    }

    /// @brief Convenience constructor using Rivet enums for most common jet algs (including some plugins).
    ///
    /// For the built-in algs, E-scheme recombination is used. For full control
    /// of FastJet built-in jet algs, use the constructors from native-args or a
    /// plugin pointer.
    ///
    /// @warning Provided area definition pointer must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             Algo alg, double rparameter,
             JetAlg::Muons usemuons=JetAlg::Muons::ALL,
             JetAlg::Invisibles useinvis=JetAlg::Invisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr,
             double seed_threshold=1.0)
      : JetAlg(fsp, usemuons, useinvis)
    {
      _initBase();
      _initJdef(alg, rparameter, seed_threshold);
    }


    // /// Same thing as above, but without an FS (for when we want to pass the particles directly to the calc method)
    // /// @todo Does this work properly, without internal HeavyQuarks etc.?
    // FastJets(Algo alg, double rparameter, double seed_threshold=1.0) { _initJdef(alg, rparameter, seed_threshold); }
    // /// Same thing as above, but without an FS (for when we want to pass the particles directly to the calc method)
    // /// @todo Does this work properly, without internal HeavyQuarks etc.?
    // FastJets(fastjet::JetAlgorithm type, fastjet::RecombinationScheme recom, double rparameter) { _initJdef(type, recom, rparameter); }
    // /// Same thing as above, but without an FS (for when we want to pass the particles directly to the calc method)
    // /// @todo Does this work properly, without internal HeavyQuarks etc.?
    // FastJets(fastjet::JetDefinition::Plugin* plugin) : _jdef(plugin), _plugin(plugin) {
    //   // _plugin.reset(plugin);
    //   // _jdef = fastjet::JetDefinition(plugin);
    // }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(FastJets);

    /// @}


    /// @name Static helper functions for FastJet interaction, with tagging
    /// @{

    /// Make PseudoJets for input to a ClusterSequence, with user_index codes for constituent- and tag-particle linking
    static PseudoJets mkClusterInputs(const Particles& fsparticles, const Particles& tagparticles=Particles());
    /// Make a Rivet Jet from a PseudoJet holding a user_index code for lookup of Rivet fsparticle or tagparticle links
    static Jet mkJet(const PseudoJet& pj, const Particles& fsparticles, const Particles& tagparticles=Particles());
    /// Convert a whole list of PseudoJets to a list of Jets, with mkJet-style unpacking
    static Jets mkJets(const PseudoJets& pjs, const Particles& fsparticles, const Particles& tagparticles=Particles());

    /// @}


    /// Reset the projection. Jet def, etc. are unchanged.
    void reset();


    /// @name Jet-area calculations
    /// @{

    /// @brief Use provided jet area definition
    ///
    /// @warning The provided pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    /// @note Provide an adef null pointer to re-disable jet area calculation
    void useJetArea(fastjet::AreaDefinition* adef) {
      _adef.reset(adef);
    }

    /// Don't calculate a jet area
    void clearJetArea() {
      _adef.reset();
    }

    /// @}


    /// @name Jet grooming
    /// @{

    /// @brief Add a grooming transformer (base class of fastjet::Filter, etc.)
    ///
    /// @warning The provided pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    /// @note Provide an adef null pointer to re-disable jet area calculation
    void addTrf(fastjet::Transformer* trf) {
      _trfs.push_back(shared_ptr<fastjet::Transformer>(trf));
    }

    /// @brief Add a list of grooming transformers
    ///
    /// @warning The provided pointers must be heap-allocated: they will be stored/deleted via a shared_ptr.
    /// @note Provide an adef null pointer to re-disable jet area calculation
    template<typename TRFS, typename TRF=typename TRFS::value_type>
    typename std::enable_if<Derefable<TRF>::value, void>::type
    addTrfs(const TRFS& trfs) {
      for (auto& trf : trfs) addTrf(trf);
    }

    /// Don't apply any jet transformers
    void clearTrfs() {
      _trfs.clear();
    }

    /// @brief Trim (filter) a jet, keeping tag and constituent info in the resulting jet
    ///
    /// @deprecated Use the built-in transformers system, e.g. addTrf(), instead
    Jet trimJet(const Jet& input, const fastjet::Filter& trimmer) const;

    /// @}


    /// @name Access to the jets
    //@{

    /// Get the jets (unordered) with pT > ptmin.
    Jets _jets() const;

    /// Get the pseudo jets (unordered).
    /// @deprecated Use pseudojets
    PseudoJets pseudoJets(double ptmin=0.0) const;
    /// Alias
    PseudoJets pseudojets(double ptmin=0.0) const { return pseudoJets(ptmin); }

    /// Get the pseudo jets, ordered by \f$ p_T \f$.
    /// @deprecated Use pseudojetsbyPt
    PseudoJets pseudoJetsByPt(double ptmin=0.0) const {
      return sorted_by_pt(pseudoJets(ptmin));
    }
    /// Alias
    PseudoJets pseudojetsByPt(double ptmin=0.0) const { return pseudoJetsByPt(ptmin); }

    /// Get the pseudo jets, ordered by \f$ E \f$.
    /// @deprecated Use pseudojetsByE
    PseudoJets pseudoJetsByE(double ptmin=0.0) const {
      return sorted_by_E(pseudoJets(ptmin));
    }
    /// Alias
    PseudoJets pseudojetsByE(double ptmin=0.0) const { return pseudoJetsByE(ptmin); }

    /// Get the pseudo jets, ordered by rapidity.
    /// @deprecated Use pseudojetsByRapidity
    PseudoJets pseudoJetsByRapidity(double ptmin=0.0) const {
      return sorted_by_rapidity(pseudoJets(ptmin));
    }
    /// Alias
    PseudoJets pseudojetsByRapidity(double ptmin=0.0) const { return pseudoJetsByRapidity(ptmin); }

    //@}


    /// @name Access to the FastJet clustering objects such as jet def, area def, and cluster
    //@{

    /// Return the cluster sequence.
    /// @todo Care needed re. const shared_ptr<T> vs. shared_ptr<const T>
    const shared_ptr<fastjet::ClusterSequence> clusterSeq() const {
      return _cseq;
    }

    /// Return the area-enabled cluster sequence (if an area defn exists, otherwise returns a null ptr).
    /// @todo Care needed re. const shared_ptr<T> vs. shared_ptr<const T>
    const shared_ptr<fastjet::ClusterSequenceArea> clusterSeqArea() const {
      return areaDef() ? dynamic_pointer_cast<fastjet::ClusterSequenceArea>(_cseq) : nullptr;
    }

    /// Return the jet definition.
    const fastjet::JetDefinition& jetDef() const {
      return _jdef;
    }

    /// @brief Return the area definition.
    ///
    /// @warning May be null!
    /// @todo Care needed re. const shared_ptr<T> vs. shared_ptr<const T>
    const shared_ptr<fastjet::AreaDefinition> areaDef() const {
      return _adef;
    }

    //@}


  private:

    /// Shared utility functions to implement constructor behaviour
    void _initBase();
    void _initJdef(Algo alg, double rparameter, double seed_threshold);

  protected:

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;

  public:

    /// Do the calculation locally (no caching).
    void calc(const Particles& fsparticles, const Particles& tagparticles=Particles());


  private:

    /// Jet definition
    fastjet::JetDefinition _jdef;

    /// Pointer to user-handled area definition
    std::shared_ptr<fastjet::AreaDefinition> _adef;

    /// Cluster sequence
    std::shared_ptr<fastjet::ClusterSequence> _cseq;

    /// FastJet external plugin
    std::shared_ptr<fastjet::JetDefinition::Plugin> _plugin;

    /// List of jet groomers to be applied
    std::vector< std::shared_ptr<fastjet::Transformer> > _trfs;

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable std::map<int, vector<double> > _yscales;

    /// Particles used for constituent and tag lookup
    Particles _fsparticles, _tagparticles;

  };

}

#endif

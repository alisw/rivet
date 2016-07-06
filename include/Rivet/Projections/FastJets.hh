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
//#include "fastjet/PxConePlugin.hh"

namespace Rivet {


  /// Project out jets found using the FastJet package jet algorithms.
  class FastJets : public JetAlg {
  public:

    /// Wrapper enum for selected Fastjet jet algorithms.
    enum JetAlgName { KT, CAM, SISCONE, ANTIKT,
                      PXCONE,
                      ATLASCONE, CMSCONE,
                      CDFJETCLU, CDFMIDPOINT, D0ILCONE,
                      JADE, DURHAM, TRACKJET };


    /// @name Constructors etc.
    //@{

    /// "Wrapped" argument constructor using Rivet enums for most common
    /// jet alg choices (including some plugins). For the built-in algs,
    /// E-scheme recombination is used. For full control of
    /// FastJet built-in jet algs, use the native arg constructor.
    FastJets(const FinalState& fsp, JetAlgName alg,
             double rparameter, double seed_threshold=1.0)
      : JetAlg(fsp), _adef(0) { _init1(alg, rparameter, seed_threshold); }

    /// Native argument constructor, using FastJet alg/scheme enums.
    FastJets(const FinalState& fsp, fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter)
      : JetAlg(fsp), _adef(0) { _init2(type, recom, rparameter); }

    /// Explicitly pass in an externally-constructed plugin (must be heap-allocated, Rivet will delete)
    FastJets(const FinalState& fsp, fastjet::JetDefinition::Plugin* plugin)
      : JetAlg(fsp), _adef(0) { _init3(plugin); }
    /// Explicitly pass in an externally-constructed plugin (must be heap-allocated, Rivet will delete)
    FastJets(const FinalState& fsp, fastjet::JetDefinition::Plugin& plugin)
      : JetAlg(fsp), _adef(0) { _init3(&plugin); }


    /// Same thing as above, but without an FS (for when we want to pass the particles directly to the calc method)
    FastJets(JetAlgName alg, double rparameter, double seed_threshold=1.0)
      : _adef(0) { _init1(alg, rparameter, seed_threshold); }
    /// Same thing as above, but without an FS (for when we want to pass the particles directly to the calc method)
    FastJets(fastjet::JetAlgorithm type, fastjet::RecombinationScheme recom, double rparameter)
      : _adef(0) { _init2(type, recom, rparameter); }
    /// Same thing as above, but without an FS (for when we want to pass the particles directly to the calc method)
    FastJets(fastjet::JetDefinition::Plugin* plugin)
      : _adef(0) { _init3(plugin); }
    /// Same thing as above, but without an FS (for when we want to pass the particles directly to the calc method)
    FastJets(fastjet::JetDefinition::Plugin& plugin)
      : _adef(0) { _init3(&plugin); }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new FastJets(*this);
    }

    //@}


  public:

    /// Reset the projection. Jet def, etc. are unchanged.
    void reset();

    /// @brief Use provided jet area definition
    ///
    /// @warning @a adef is NOT copied, the user must ensure that it remains valid!
    ///
    /// Provide an adef null pointer to re-disable jet area calculation
    void useJetArea(fastjet::AreaDefinition* adef) {
      _adef = adef;
    }

    /// @name Access to the jets
    //@{

    /// Number of jets above the \f$ p_\perp \f$ cut.
    size_t numJets(double ptmin=0.0) const;

    /// Number of jets.
    size_t size() const {
      return numJets();
    }

    /// Get the jets (unordered).
    Jets _jets(double ptmin=0.0) const;

    /// Get the pseudo jets (unordered).
    PseudoJets pseudoJets(double ptmin=0.0) const;
    /// Alias
    PseudoJets pseudojets(double ptmin=0.0) const { return pseudoJets(ptmin); }

    /// Get the pseudo jets, ordered by \f$ p_T \f$.
    PseudoJets pseudoJetsByPt(double ptmin=0.0) const {
      return sorted_by_pt(pseudoJets(ptmin));
    }
    /// Alias
    PseudoJets pseudojetsByPt(double ptmin=0.0) const { return pseudoJetsByPt(ptmin); }

    /// Get the pseudo jets, ordered by \f$ E \f$.
    PseudoJets pseudoJetsByE(double ptmin=0.0) const {
      return sorted_by_E(pseudoJets(ptmin));
    }
    /// Alias
    PseudoJets pseudojetsByE(double ptmin=0.0) const { return pseudoJetsByE(ptmin); }

    /// Get the pseudo jets, ordered by rapidity.
    PseudoJets pseudoJetsByRapidity(double ptmin=0.0) const {
      return sorted_by_rapidity(pseudoJets(ptmin));
    }
    /// Alias
    PseudoJets pseudojetsByRapidity(double ptmin=0.0) const { return pseudoJetsByRapidity(ptmin); }

    //@}


    /// @name Access to the cluster sequence and splitting info
    //@{

    /// Return the cluster sequence.
    const fastjet::ClusterSequence* clusterSeq() const {
      return _cseq.get();
    }

    /// Return the area-enabled cluster sequence (if an area defn exists, otherwise returns a null ptr).
    const fastjet::ClusterSequenceArea* clusterSeqArea() const {
      if (areaDef() == NULL) return (fastjet::ClusterSequenceArea*) 0;
      return dynamic_cast<fastjet::ClusterSequenceArea*>(_cseq.get());
    }

    /// Return the jet definition.
    const fastjet::JetDefinition& jetDef() const {
      return _jdef;
    }

    /// Return the area definition.. May be null.
    const fastjet::AreaDefinition* areaDef() const {
      return _adef;
    }


    /// Get the subjet splitting variables for the given jet.
    vector<double> ySubJet(const fastjet::PseudoJet& jet) const;

    /// @brief Split a jet a la PRL100,242001(2008).
    ///
    /// Based on code from G.Salam, A.Davison.
    fastjet::PseudoJet splitJet(fastjet::PseudoJet jet, double& last_R) const;

    /// @brief Filter a jet a la PRL100,242001(2008).
    ///
    /// Based on code from G.Salam, A.Davison.
    fastjet::PseudoJet filterJet(fastjet::PseudoJet jet, double& stingy_R, const double def_R) const;

    //@}


  private:

    /// Shared utility functions to implement constructor behaviour
    void _init1(JetAlgName alg, double rparameter, double seed_threshold);
    void _init2(fastjet::JetAlgorithm type, fastjet::RecombinationScheme recom, double rparameter);
    void _init3(fastjet::JetDefinition::Plugin* plugin);

  protected:

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;

  public:

    /// Do the calculation locally (no caching).
    void calc(const Particles& fsparticles, const Particles& tagparticles=Particles());


  private:

    /// Jet definition
    fastjet::JetDefinition _jdef;

    /// Pointer to user-handled area definition
    fastjet::AreaDefinition* _adef;

    /// Cluster sequence
    shared_ptr<fastjet::ClusterSequence> _cseq;

    /// FastJet external plugin
    shared_ptr<fastjet::JetDefinition::Plugin> _plugin;

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable map<int, vector<double> > _yscales;

    /// set of particles sorted by their PT2
    //set<Particle, ParticleBase::byPTAscending> _particles;
    map<int, Particle> _particles;

  };

}

#endif

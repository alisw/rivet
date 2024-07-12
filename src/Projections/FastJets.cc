// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  void FastJets::_initBase() {
    setName("FastJets");
    declare(HeavyHadrons(), "HFHadrons");
    declare(TauFinder(TauFinder::DecayMode::HADRONIC), "Taus");

    // Print/hide FJ banner
    std::cout.setstate(std::ios_base::badbit);
    fastjet::ClusterSequence::print_banner();
    std::cout.clear();
  }


  void FastJets::_initJdef(Algo alg, double rparameter, double seed_threshold) {
    MSG_DEBUG("JetAlg = " << static_cast<int>(alg));
    MSG_DEBUG("R parameter = " << rparameter);
    MSG_DEBUG("Seed threshold = " << seed_threshold);
    if (alg == KT) {
      _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == ANTIKT) {
      _jdef = fastjet::JetDefinition(fastjet::antikt_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == CAM) {
      _jdef = fastjet::JetDefinition(fastjet::cambridge_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == DURHAM) {
      _jdef = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::E_scheme);
    } else if (alg == GENKTEE) {
      _jdef = fastjet::JetDefinition(fastjet::ee_genkt_algorithm, rparameter, -1);
    } else if (alg == KTET) {
      _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, rparameter, fastjet::Et_scheme);
    } else if (alg == ANTIKTET) {
      _jdef = fastjet::JetDefinition(fastjet::antikt_algorithm, rparameter, fastjet::Et_scheme);

    } else {
      // Plugins:
      if (alg == SISCONE) {
        const double OVERLAP_THRESHOLD = 0.75;
        _plugin.reset(new fastjet::SISConePlugin(rparameter, OVERLAP_THRESHOLD));
      } else if (alg == PXCONE) {
        string msg = "Using own c++ version of PxCone, since FastJet doesn't install it by default. ";
        msg += "Please notify the Rivet authors if this behaviour should be changed.";
        MSG_WARNING(msg);
        _plugin.reset(new Rivet::PxConePlugin(rparameter));
      } else if (alg == ATLASCONE) {
        const double OVERLAP_THRESHOLD = 0.5;
        _plugin.reset(new fastjet::ATLASConePlugin(rparameter, seed_threshold, OVERLAP_THRESHOLD));
      } else if (alg == CMSCONE) {
        _plugin.reset(new fastjet::CMSIterativeConePlugin(rparameter, seed_threshold));
      } else if (alg == CDFJETCLU) {
        const double OVERLAP_THRESHOLD = 0.75;
        _plugin.reset(new fastjet::CDFJetCluPlugin(rparameter, OVERLAP_THRESHOLD, seed_threshold));
      } else if (alg == CDFMIDPOINT) {
        const double OVERLAP_THRESHOLD = 0.5;
        _plugin.reset(new fastjet::CDFMidPointPlugin(rparameter, OVERLAP_THRESHOLD, seed_threshold));
      } else if (alg == D0ILCONE) {
        const double min_jet_Et = 6.0;
        _plugin.reset(new fastjet::D0RunIIConePlugin(rparameter, min_jet_Et));
      } else if (alg == JADE) {
        _plugin.reset(new fastjet::JadePlugin());
      } else if (alg == TRACKJET) {
        _plugin.reset(new fastjet::TrackJetPlugin(rparameter));
      }
      _jdef = fastjet::JetDefinition(_plugin.get());
    }
  }


  CmpState FastJets::compare(const Projection& p) const {
    const FastJets& other = dynamic_cast<const FastJets&>(p);
    CmpState rtn =
      cmp(_useMuons, other._useMuons) ||
      cmp(_useInvisibles, other._useInvisibles) ||
      mkNamedPCmp(other, "FS") ||
      cmp(_jdef.jet_algorithm(), other._jdef.jet_algorithm()) ||
      cmp(_jdef.recombination_scheme(), other._jdef.recombination_scheme()) ||
      cmp(_jdef.plugin(), other._jdef.plugin()) ||
      cmp(_jdef.R(), other._jdef.R()) ||
      cmp(_adef, other._adef);
    if (rtn != CmpState::EQ) return rtn; //< shortcut transformer comparison if aleady different

    // Compare the transformers list
    if (_trfs.empty() && other._trfs.empty()) return CmpState::EQ;
    /// @todo Improve fastjet::Transformer to add a virtual operator==, and use all()
    // if (_trfs.size() != other._trfs.size()) return CmpState::NEQ;
    // for (size_t it = 0; it < _trfs.size(); ++it) {

    return CmpState::NEQ;
  }


  // STATIC
  PseudoJets FastJets::mkClusterInputs(const Particles& fsparticles, const Particles& tagparticles) {
    PseudoJets pjs;
    /// @todo Use FastJet3's UserInfo system to store Particle pointers directly?

    // Store 4 vector data about each particle into FastJet's PseudoJets
    for (size_t i = 0; i < fsparticles.size(); ++i) {
      fastjet::PseudoJet pj = fsparticles[i];
      pj.set_user_index(i+1);
      pjs.push_back(pj);
    }
    // And the same for ghost tagging particles (with negative user indices)
    for (size_t i = 0; i < tagparticles.size(); ++i) {
      fastjet::PseudoJet pj = tagparticles[i];
      pj *= 1e-20; ///< Ghostify the momentum
      pj.set_user_index(-i-1);
      pjs.push_back(pj);
    }

    return pjs;
  }


  // STATIC
  Jet FastJets::mkJet(const PseudoJet& pj, const Particles& fsparticles, const Particles& tagparticles) {
    const PseudoJets pjconstituents = pj.constituents();

    Particles constituents, tags;
    constituents.reserve(pjconstituents.size());

    for (const fastjet::PseudoJet& pjc : pjconstituents) {
      // Pure ghosts don't have corresponding particles
      if (pjc.has_area() && pjc.is_pure_ghost()) continue;
      // Default user index = 0 doesn't give valid particle lookup
      if (pjc.user_index() == 0) continue;
      // Split by index sign into constituent & tag lookup
      if (pjc.user_index() > 0) {
        // Find constituents if index > 0
        const size_t i = pjc.user_index() - 1;
        if (i >= fsparticles.size()) throw RangeError("FS particle lookup failed in jet construction");
        constituents.push_back(fsparticles.at(i));
      } else if (!tagparticles.empty()) {
        // Find tags if index < 0
        const size_t i = abs(pjc.user_index()) - 1;
        if (i >= tagparticles.size()) throw RangeError("Tag particle lookup failed in jet construction");
        tags.push_back(tagparticles.at(i));
      }
    }

    return Jet(pj, constituents, tags);
  }


  // STATIC
  Jets FastJets::mkJets(const PseudoJets& pjs, const Particles& fsparticles, const Particles& tagparticles) {
    Jets rtn; rtn.reserve(pjs.size());
    for (const PseudoJet& pj : pjs) {
      rtn.push_back(FastJets::mkJet(pj, fsparticles, tagparticles));
    }
    return rtn;
  }


  void FastJets::project(const Event& e) {
    // Assemble final state particles
    const string fskey = (_useInvisibles == JetAlg::Invisibles::NONE) ? "VFS" : "FS";
    Particles fsparticles = applyProjection<FinalState>(e, fskey).particles();
    // Remove prompt invisibles if needed (already done by VFS if using NO_INVISIBLES)
    if (_useInvisibles == JetAlg::Invisibles::DECAY) {
      ifilter_discard(fsparticles, [](const Particle& p) { return !p.isVisible() && p.isPrompt(); });
    }
    // Remove prompt/all muons if needed
    if (_useMuons == JetAlg::Muons::DECAY) {
      ifilter_discard(fsparticles, [](const Particle& p) { return isMuon(p) && p.isPrompt(); });
    } else if (_useMuons == JetAlg::Muons::NONE) {
      ifilter_discard(fsparticles, isMuon);
    }

    // Tagging particles
    const Particles chadrons = applyProjection<HeavyHadrons>(e, "HFHadrons").cHadrons();
    const Particles bhadrons = applyProjection<HeavyHadrons>(e, "HFHadrons").bHadrons();
    const Particles taus = applyProjection<FinalState>(e, "Taus").particles();

    // Run the calculation
    calc(fsparticles, chadrons+bhadrons+taus);
  }


  void FastJets::calc(const Particles& fsparticles, const Particles& tagparticles) {
    MSG_DEBUG("Finding jets from " << fsparticles.size() << " input particles + " << tagparticles.size() << " tagging particles");
    _fsparticles = fsparticles;
    _tagparticles = tagparticles;

    // Make pseudojets, with mapping info to Rivet FS and tag particles
    PseudoJets pjs = mkClusterInputs(_fsparticles, _tagparticles);

    // Run either basic or area-calculating cluster sequence as reqd.
    /// @todo How can we make sure these persist if we call calc multiple times?
    if (_adef) {
      _cseq.reset(new fastjet::ClusterSequenceArea(pjs, _jdef, *_adef));
    } else {
      _cseq.reset(new fastjet::ClusterSequence(pjs, _jdef));
    }

    MSG_DEBUG("ClusterSequence constructed; Njets_tot = "
              << _cseq->inclusive_jets().size() << ", Njets(pT > 10 GeV) = "
              << _cseq->inclusive_jets(10*GeV).size());
  }


  void FastJets::reset() {
    _yscales.clear();
    _fsparticles.clear();
    _tagparticles.clear();
    /// @todo _cseq = fastjet::ClusterSequence();
  }


  Jets FastJets::_jets() const {
    /// @todo Cache?
    return mkJets(pseudojets(), _fsparticles, _tagparticles);
  }


  /// @todo "Automate" trimming as part of project() with pre-registered Filters
  Jet FastJets::trimJet(const Jet& input, const fastjet::Filter& trimmer) const {
    if (input.pseudojet().associated_cluster_sequence() != clusterSeq().get())
      throw Error("To trim a Rivet::Jet, its associated PseudoJet must have come from this FastJets' ClusterSequence");
    PseudoJet pj = trimmer(input);
    return mkJet(pj, _fsparticles, _tagparticles);
  }


  PseudoJets FastJets::pseudoJets(double ptmin) const {
    // Get the base set of pseudo-jets
    PseudoJets rtn = clusterSeq() ? clusterSeq()->inclusive_jets(ptmin) : PseudoJets();

    // Run the jet groomers on each jet
    for (PseudoJet& pj : rtn) {
      for (auto& t : _trfs) {
        pj = t->result(pj);
      }
    }

    return rtn;
  }


}

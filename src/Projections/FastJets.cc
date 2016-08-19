// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  void FastJets::_initBase() {
    setName("FastJets");
    addProjection(HeavyHadrons(), "HFHadrons");
    addProjection(TauFinder(TauFinder::HADRONIC), "Taus");
  }


  // void FastJets::_init1(JetAlgName alg, double rparameter, double seed_threshold) {
  //   _initBase();


  void FastJets::_initJdef(JetAlgName alg, double rparameter, double seed_threshold) {
    MSG_DEBUG("JetAlg = " << alg);
    MSG_DEBUG("R parameter = " << rparameter);
    MSG_DEBUG("Seed threshold = " << seed_threshold);
    if (alg == KT) {
      _jdef = fastjet::JetDefinition(fastjet::kt_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == CAM) {
      _jdef = fastjet::JetDefinition(fastjet::cambridge_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == ANTIKT) {
      _jdef = fastjet::JetDefinition(fastjet::antikt_algorithm, rparameter, fastjet::E_scheme);
    } else if (alg == DURHAM) {
      _jdef = fastjet::JetDefinition(fastjet::ee_kt_algorithm, fastjet::E_scheme);
    } else {
      // Plugins:
      if (alg == SISCONE) {
        const double OVERLAP_THRESHOLD = 0.75;
        _plugin.reset(new fastjet::SISConePlugin(rparameter, OVERLAP_THRESHOLD));
      // } else if (alg == PXCONE) {
      //   string msg = "PxCone currently not supported, since FastJet doesn't install it by default. ";
      //   msg += "Please notify the Rivet authors if this behaviour should be changed.";
      //   throw Error(msg);
      //  _plugin.reset(new fastjet::PxConePlugin(rparameter));
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


  // void FastJets::_init2(fastjet::JetAlgorithm type,
  //                       fastjet::RecombinationScheme recom, double rparameter)
  //   : FastJets
  // {

  //   _initBase();
  //   _jdef = fastjet::JetDefinition(type, rparameter, recom);
  // }


  // void FastJets::_init3(const fastjet::JetDefinition& jdef) {
  //   _initBase();
  //   _jdef = jdef;
  // }


  // void FastJets::_init4(fastjet::JetDefinition::Plugin* plugin) {
  //   _initBase();
  //   _plugin.reset(plugin);
  //   _jdef = fastjet::JetDefinition(_plugin.get());
  // }



  int FastJets::compare(const Projection& p) const {
    const FastJets& other = dynamic_cast<const FastJets&>(p);
    return \
      cmp(_useMuons, other._useMuons) ||
      cmp(_useInvisibles, other._useInvisibles) ||
      mkNamedPCmp(other, "FS") ||
      cmp(_jdef.jet_algorithm(), other._jdef.jet_algorithm()) ||
      cmp(_jdef.recombination_scheme(), other._jdef.recombination_scheme()) ||
      cmp(_jdef.plugin(), other._jdef.plugin()) ||
      cmp(_jdef.R(), other._jdef.R()) ||
      cmp(_adef, other._adef);
  }


  namespace {
    /// @todo Replace with C++11 lambdas
    bool isPromptInvisible(const Particle& p) { return !(p.isVisible() || p.fromDecay()); }
    // bool isMuon(const Particle& p) { return p.abspid() == PID::MUON; }
    bool isPromptMuon(const Particle& p) { return isMuon(p) && !p.fromDecay(); }
  }


  void FastJets::project(const Event& e) {
    // Assemble final state particles
    const string fskey = (_useInvisibles == JetAlg::NO_INVISIBLES) ? "VFS" : "FS";
    Particles fsparticles = applyProjection<FinalState>(e, fskey).particles();
    // Remove prompt invisibles if needed (already done by VFS if using NO_INVISIBLES)
    if (_useInvisibles == JetAlg::DECAY_INVISIBLES)
      fsparticles.erase( std::remove_if(fsparticles.begin(), fsparticles.end(), isPromptInvisible), fsparticles.end() );
    // Remove prompt/all muons if needed
    if (_useMuons == JetAlg::DECAY_MUONS)
      fsparticles.erase( std::remove_if(fsparticles.begin(), fsparticles.end(), isPromptMuon), fsparticles.end() );
    else if (_useMuons == JetAlg::NO_MUONS)
      fsparticles.erase( std::remove_if(fsparticles.begin(), fsparticles.end(), isMuon), fsparticles.end() );

    // Tagging particles
    /// @todo Allow the user to specify tag particle kinematic thresholds
    const Particles chadrons = applyProjection<HeavyHadrons>(e, "HFHadrons").cHadrons();
    const Particles bhadrons = applyProjection<HeavyHadrons>(e, "HFHadrons").bHadrons();
    const Particles taus = applyProjection<FinalState>(e, "Taus").particles();
    calc(fsparticles, chadrons+bhadrons+taus);
  }


  void FastJets::calc(const Particles& fsparticles, const Particles& tagparticles) {
    _particles.clear();
    vector<fastjet::PseudoJet> pjs;

    MSG_DEBUG("Finding jets from " << fsparticles.size() << " input particles");

    /// @todo Use FastJet3's UserInfo system

    // Store 4 vector data about each particle into FastJet's PseudoJets
    int counter = 1;
    for (const Particle& p : fsparticles) {
      fastjet::PseudoJet pj = p;
      pj.set_user_index(counter);
      pjs.push_back(pj);
      _particles[counter] = p;
      counter += 1;
    }
    // And the same for ghost tagging particles (with negative user indices)
    counter = 1;
    for (const Particle& p : tagparticles) {
      fastjet::PseudoJet pj = p;
      pj *= 1e-20; ///< Ghostify the momentum
      pj.set_user_index(-counter);
      pjs.push_back(pj);
      _particles[-counter] = p;
      counter += 1;
    }

    // Choose cseq as basic or area-calculating
    if (_adef) {
      _cseq.reset(new fastjet::ClusterSequenceArea(pjs, _jdef, *_adef));
    } else {
      _cseq.reset(new fastjet::ClusterSequence(pjs, _jdef));
    }
    MSG_DEBUG("FastJet ClusterSequence constructed; Njets_tot = "
              << _cseq->inclusive_jets().size() << ", Njets_10 = "
              << _cseq->inclusive_jets(10*GeV).size()); //< only inefficient in debug mode
  }


  void FastJets::reset() {
    _yscales.clear();
    _particles.clear();
    /// @todo _cseq = fastjet::ClusterSequence();
  }


  Jets FastJets::_jets() const {
    Jets rtn; rtn.reserve(pseudojets().size());
    foreach (const fastjet::PseudoJet& pj, pseudojets()) {
      rtn.push_back(_mkJet(pj));
    }
    /// @todo Cache?
    return rtn;
  }


  Jet FastJets::trimJet(const Jet &input, const fastjet::Filter &trimmer)const{
    assert(input.pseudojet().associated_cluster_sequence() == clusterSeq().get());
    PseudoJet pj = trimmer(input);
    return _mkJet(pj);
  }


  PseudoJets FastJets::pseudoJets(double ptmin) const {
    return clusterSeq() ? clusterSeq()->inclusive_jets(ptmin) : PseudoJets();
  }


  Jet FastJets::_mkJet(const PseudoJet &pj)const{
    assert(clusterSeq());

    // Take the constituents from the cluster sequence, unless the jet was not
    // associated with the cluster sequence (could be the case for trimmed jets)
    const PseudoJets parts = (pj.associated_cluster_sequence() == clusterSeq().get())
      ? clusterSeq()->constituents(pj) : pj.constituents();

    vector<Particle> constituents, tags;
    constituents.reserve(parts.size());
    for (const fastjet::PseudoJet& p : parts) {
      map<int, Particle>::const_iterator found = _particles.find(p.user_index());
      // assert(found != _particles.end());
      if (found == _particles.end() && p.is_pure_ghost()) continue; //< Pure FJ ghosts are ok
      assert(found != _particles.end()); //< Anything else must be known
      assert(found->first != 0); //< All mapping IDs are pos-def (particles) or neg-def (tags)
      if (found->first > 0) constituents.push_back(found->second);
      else if (found->first < 0) tags.push_back(found->second);
    }

    return Jet(pj, constituents, tags);
  }


}

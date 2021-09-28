// -*- C++ -*-
#ifndef RIVET_EventMixingFinalState_HH
#define RIVET_EventMixingFinalState_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/Random.hh"
#include <deque>
#include <algorithm>

namespace Rivet {


  // @brief Projects out an event mixed of several events, given
  // a mixing observable (eg. number of final state particles),
  // defining what should qualify as a mixable event.
  // Binning in the mixing observable is defined in the constructor,
  // as is the number of events one wants to mix with.
  // The method calculateMixingObs() must can be overloaded
  // in derived classes, to provide the definition of the mixing observable,
  // on the provided projection, eg. centrality or something more elaborate.
  //
  // The implementation can correcly handle mixing of weighted events, but
  // not multi-weighted events. Nor does the implementation attemt to handle
  // calculation of one (or more) event weights for the mixed events. For most
  // common use cases, like calculating a background sample, this is sufficient.
  // If a more elaborate use case ever turns up, this must be reevaluated.
  //
  //
  // @author Christian Bierlich <christian.bierlich@thep.lu.se>

  // Weighted random shuffle, similar to std::random_shuffle, which
  // allows the passing of a weight for each element to be shuffled.
  template <class RandomAccessIterator,
            class WeightIterator, class RandomNumberGenerator>
  void weighted_shuffle(RandomAccessIterator first, RandomAccessIterator last,
                        WeightIterator fw, WeightIterator lw, RandomNumberGenerator& g) {
    while(first != last && fw != lw) {
      std::discrete_distribution<int> weightDist(fw, lw);
      int i = weightDist(g);
      if(i){
        std::iter_swap(first, next(first, i));
        std::iter_swap(fw, next(fw, i));
      }
      ++first;
      ++fw;
    }
  }
  // A MixEvent is a vector of particles with and associated weight.
  typedef pair<Particles, double> MixEvent;
  typedef map<double, std::deque<MixEvent> > MixMap;

  /// EventMixingBase is the base class for event mixing projections.
  ///
  /// Most methods are defined in this base class as they should.
  /// In order to use it, a derived class should be implemented where:
  /// - The constructor is reimplmented, giving the derived projection type
  /// from which the mixing observable is calculated. The constructor must also
  /// be declared public in the derived class.
  /// - The calculateMixingObs is implemented.
  /// To examples of such derived classes are given below,
  /// 1) EventMixingFinalState, where the mixing observable are calculated
  /// on a multiplicity of a charged final state, and:
  /// 2) EventMixingCentrality, where the mixing observable is centrality.
  ///
  class EventMixingBase : public Projection {
  protected:
    // Constructor
    EventMixingBase(const Projection & mixObsProj, const ParticleFinder& mix,
                    size_t nMixIn, double oMin, double oMax, double deltao,
                    const size_t defaultIdx) : nMix(nMixIn), unitWeights(true) {
      // The base class contructor should be called explicitly in derived classes
      // to add projections below.
      setName("EventMixingBase");
      declare(mixObsProj,"OBS");
      declare(mix,"MIX");
      MSG_WARNING("EventMixing is not fully validated. Use with caution.");

      _defaultWeightIdx = defaultIdx;
      // Set up the map for mixing events.
      for(double o = oMin; o < oMax; o+=deltao )
        mixEvents[o] = std::deque<MixEvent>();
    }

  public:

    // Test if we have enough mixing events available for projected,
    // current mixing observable.
    bool hasMixingEvents() const {
      MixMap::const_iterator mixItr = mixEvents.lower_bound(mObs);
      if(mixItr == mixEvents.end() || mixItr->second.size() < nMix + 1)
        return false;
      return true;
    }

    // Return a vector of mixing events.
    vector<MixEvent> getMixingEvents() const {
      if (!hasMixingEvents())
        return vector<MixEvent>();
      MixMap::const_iterator mixItr = mixEvents.lower_bound(mObs);
      return vector<MixEvent>(mixItr->second.begin(), mixItr->second.end() - 1);
    }

    // Return a vector of particles from the mixing events. Can
    // be overloaded in derived classes, though normally not neccesary.
    virtual const Particles particles() const {
      // Test if we have enough mixing events.
      if (!hasMixingEvents())
        return Particles();
      // Get mixing events for the current, projected mixing observable.
      MixMap::const_iterator mixItr = mixEvents.lower_bound(mObs);
      vector<MixEvent> mixEvents(mixItr->second.begin(), mixItr->second.end() - 1);
      // Make the vector of mixed particles.
      Particles mixParticles;
      vector<double> weights;
      size_t pSize = 0;
      for (size_t i = 0; i < mixEvents.size(); ++i)
        pSize+=mixEvents[i].first.size();
      mixParticles.reserve(pSize);
      weights.reserve(pSize);
      // Put the particles in the vector.
      for (size_t i = 0; i < mixEvents.size(); ++i) {
        mixParticles.insert(mixParticles.end(), mixEvents[i].first.begin(), mixEvents[i].first.end());
        vector<double> tmp(mixEvents[i].first.size(), mixEvents[i].second);
        weights.insert(weights.end(), tmp.begin(), tmp.end());
      }

      // Shuffle the particles.
      if (unitWeights) {
        // Use the thread safe random number generator.
        //auto rnd = [&] (int i) {return rng()()%i;};
        std::shuffle(mixParticles.begin(), mixParticles.end(), rng());
        return mixParticles;
      } else {
        weighted_shuffle(mixParticles.begin(), mixParticles.end(), weights.begin(), weights.end(), rng());
        Particles tmp = vector<Particle>(mixParticles.begin(), mixParticles.begin() + size_t(ceil(mixParticles.size() / 2)));
        return tmp;
      }
    }


  protected:

    // Calulate mixing observable.
    // Must be overloaded in derived classes.
    virtual void calculateMixingObs(const Projection* mProj) = 0;


    /// Perform the projection on the Event.
    void project(const Event& e) {
      const Projection* mixObsProjPtr = &applyProjection<Projection>(e, "OBS");
      calculateMixingObs(mixObsProjPtr);
      MixMap::iterator mixItr = mixEvents.lower_bound(mObs);
      if (mixItr == mixEvents.end()){
        // We are out of bounds.
        MSG_DEBUG("Mixing observable out of bounds.");
        return;
      }
      const Particles mix = applyProjection<ParticleFinder>(e, "MIX").particles();
      mixItr->second.push_back(make_pair(mix,e.weights()[_defaultWeightIdx]));
      // Assume unit weights until we see otherwise.
      if (unitWeights && e.weights()[_defaultWeightIdx] != 1.0 ) {
        unitWeights = false;
        nMix *= 2;
      }
      if (mixItr->second.size() > nMix + 1)
        mixItr->second.pop_front();
    }


    /// Compare with other projections
    CmpState compare(const Projection& p) const {
      return mkNamedPCmp(p,"OBS");
    }


    /// The mixing observable of the current event.
    double mObs;


  private:

    /// The number of event to mix with.
    size_t nMix;

    /// The event map.
    MixMap mixEvents;

    /// Using unit weights or not.
    bool unitWeights;

    size_t _defaultWeightIdx;

  };


  // EventMixingFinalState has multiplicity as the mixing observable
  class EventMixingFinalState : public EventMixingBase {
  public:
    EventMixingFinalState(const ParticleFinder & mixObsProj,
                          const ParticleFinder& mix, size_t nMixIn, double oMin, double oMax,
                          double deltao, const size_t defaultIdx) :
      EventMixingBase(mixObsProj, mix, nMixIn, oMin, oMax, deltao, defaultIdx) {
      setName("EventMixingFinalState");
    }

    DEFAULT_RIVET_PROJ_CLONE(EventMixingFinalState);

  protected:

    // Calculate mixing observable
    virtual void calculateMixingObs(const Projection* mProj) {
      mObs = ((ParticleFinder*) mProj)->particles().size();
    }

  };


  // EventMixingCentrality has centrality as the mixing observable
  class EventMixingCentrality : public EventMixingBase {
  public:

    EventMixingCentrality(const CentralityProjection& mixObsProj,
                          const ParticleFinder& mix, size_t nMixIn, double oMin, double oMax,
                          double deltao, const size_t defaultIdx) :
      EventMixingBase(mixObsProj, mix, nMixIn, oMin, oMax, deltao, defaultIdx) {
      setName("EventMixingCentrality");
    }

    DEFAULT_RIVET_PROJ_CLONE(EventMixingCentrality);

  protected:

    virtual void calculateMixingObs(const Projection* mProj) {
      mObs = ((CentralityProjection*) mProj)->operator()();
    }

  };


}

#endif

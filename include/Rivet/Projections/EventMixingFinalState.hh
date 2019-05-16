// -*- C++ -*-
#ifndef RIVET_EventMixingFinalState_HH
#define RIVET_EventMixingFinalState_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include <deque>
#include <algorithm>


namespace Rivet {

  // @brief Projects out an event mixed of several events, given
  // a mixing observable (eg. number of final state particles),
  // defining what should qualify as a mixable event.
  // Binning in the mixing observable is defined in the constructor,
  // as is the number of events one wants to mix with.
  // The protected method calculateMixingObs() can be overloaded
  // in derived classes, to define other mixing observables, eg. 
  // centrality or something even more elaborate.
  //
  // !!!DISCLAIMER!!! The projection makes no attempt at correct handling
  // of event weights - ie. what one should use as event weight for several
  // mixed events. Proceed with caution if you do not use an MC with
  // unit weights.
  //
  // @author Christian Bierlich <christian.bierlich@thep.lu.se>

  typedef map<double, deque<Particles> > MixMap;	
  class EventMixingFinalState : public Projection {
  public:
    // Constructor
    EventMixingFinalState(const ParticleFinder& fsp, const ParticleFinder& mix, size_t nMixIn,
	      double oMin, double oMax, double deltao ) : nMix(nMixIn){
    	setName("EventMixingFinalState");
	addProjection(fsp,"FS");
	addProjection(mix,"MIX");
	MSG_WARNING("EventMixingFinalState is a naive implementation, not currently " <<
		    "validated. Use with caution.");

	// Set up the map for mixing events
	for(double o = oMin; o < oMax; o+=deltao )
		mixEvents[o] = deque<Particles>();
	
    }
    // Clone on the heap
    DEFAULT_RIVET_PROJ_CLONE(EventMixingFinalState);


    // Return a vector of mixing events.
    vector<Particles> getMixingEvents() const {
	MixMap::const_iterator mixItr = mixEvents.lower_bound(mObs);
	if(mixItr == mixEvents.end() || mixItr->second.size() < nMix + 1)
		return vector<Particles>();
	return vector<Particles>(mixItr->second.begin(), mixItr->second.end() - 1);
    }
 
  protected:

    // Calulate mixing observable.
    // Can be overloaded to define more elaborate observables.
    virtual void calculateMixingObs(const Particles& parts){
    	mObs = parts.size();
    }
    
    /// Perform the projection on the Event
    virtual void project(const Event& e){
	const Particles parts = applyProjection<ParticleFinder>(e, "FS").particles();

	calculateMixingObs(parts);

	MixMap::iterator mixItr = mixEvents.lower_bound(mObs);
	if(mixItr == mixEvents.end()){
		// We are out of bounds.
		MSG_DEBUG("Mixing observable out of bounds.");
		return;
	}
	const Particles mix = applyProjection<ParticleFinder>(e, "MIX").particles();
	
	mixItr->second.push_back(mix);
	if(mixItr->second.size() > nMix + 1)
		mixItr->second.pop_front();
    }

    /// Compare with other projections
    virtual int compare(const Projection& p) const {
	return mkNamedPCmp(p,"FS");
    }
    
  protected:
    // The number of event to mix with
    size_t nMix;
    // The mixing observable of the current event
    double mObs;
    // The event map;
    MixMap mixEvents;
  };
}
#endif

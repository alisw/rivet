## Projections

Projections are arguably the most important objects in Rivet - they're certainly the part that does most of the work. A projection is an object which calculates a property or properties from an Event. Such properties can be a single number (e.g. multiplicity in a certain phase space), a single complex object (e.g. sphericity tensor) or a list of such objects (e.g. a set of boosted jets or charged particles).

In case you're confused by the name, "projection" was chosen since their task is to project out quantities of interest from the complexities of whole events - they don't strictly obey the definition of an algebraic projection operator!


### A few details about projections

All Rivet projections inherit from the abstract `Projection` class, which defines the interface elements supported by all projections. You don't need to worry about this too much - it basically forces you to supply `project()` and `compare()` methods which can be used polymorphically by Rivet's internal mechanisms.

The projections are used by Analysis objects. Internally, projections are compared via the `Projection::compare(...)` method which allows duplicate instantiations of the same projection to be avoided - as a result, the second (and so on) time a given projection is called for a particular event, it will simply return its cached value rather than repeating the computation.


### Specific projections 
Here are some examples of projections available in the current release of Rivet:

 * Beam - obtain the beam from an event
 * FinalState - get the final state particles
 * FinalStateHCM - get the final state particles, shifted to the HCM frame
 * VetoedFinalState  - get the final state particles, minus selected particle types
 * FastJets - jet algorithms accessed via FastJet
 * MissingMomentum - what it says on the tin
 * Multiplicity - count the FS particles of various types
 * Sphericity - the sphericity tensor, from which the sphericity, oblateness and planarity can be obtained
 * ParisiTensor - linearized quadratic momentum tensor, from which the Parisi C and D variables can be obtained
 * Thrust - the thrust tensor, defining the thrust axes and scalars
 * DISKinematics - kinematics of DIS events
 * DISLepton - obtain the scattered lepton in DIS events


### How to use projections

For this example we'll use the FinalState projection in an analysis. The analysis header file must include the interface of FinalState, with a `#include "Rivet/Projections/FinalState.hh"` directive.

Create and initialize projections in the analysis `init` method:
```
void init() {
  FinalState fs(-1.0, 1.0, 0.5);
  addProjection(fs, "FS");
  addProjection(Multiplicity(fs), "Mult");
  addProjection(Thrust(fs), "Thrust");
}
```

Then, when using the projection, you'll probably want to have something like this in the `analyze(...)` method:
```
void MyAnalysis::analyze(const Event& e) {
  ...
  // Project into final state
  const FinalState& fs = appl<FinalState>(e, "FS");
  ...
}
```
Note that it's good practice to make the returned projection a `const` reference -- this guarantees that the cached result will remain unchanged between repeated calls, perhaps in distinct analyses.

If using a projection inside another projection, the same applies for initialising it in the constructor of the "client" class, and the use of the projection on each event should be in the body of the `MyProjection::project(...)` method.


### Writing a projection

The best documentation when writing a projection is to look at some existing ones and use them as templates - pick a short one like `ChargedFinalState` or `Beam` since their structure is easier to see. Issues to be particularly aware of include:

 * Use without registration - it's advised that you do the actual calculation in a user-callable method called `calc`, so that the projection can be used without having to be centrally registered and attached to the event. 
 * Caching - make sure that the `project` method stores everything you need to reproduce the calculation result with minimal CPU effort. This will require some member variables explicitly put there for caching. You might want to become familiar with the `mutable` keyword, and use it **very carefully**, if the cache needs to be updated while accessing methods of a `const` projection. Make sure you reuse cached values wherever possible, by using other projections where available.
 * Comparing with other projections - the `compare` method provides a way for Rivet's internal workings to discriminate between projections which are handled polymorphically (i.e. they're all just `Projection*` as far as Rivet's type system is concerned). `compare` should return `CmpState:EQ` if two projections are effectively identical - i.e. not just the same type, but also the same configuration parameters and equivalent internal projections - and `CmpState::NEQ` if they are not identical. See some existing projections for guidance.
 * `new` and `delete` in projections - the projection and analysis destructors are pretty ineffectual things, since they **only get called at the end of the whole run**, rather than in between each event: if `new` is called during the `project` phase, `delete` had better be called in the same phase or you'll have a horrendous memory leak. Anyway, you shouldn't be using `new` unless you have a good reason, and remember: references are just as efficient as pointers and a lot safer!
As well as the logical decisions to be made in designing a projection, to integrate a projection into the build system, you will have to modify the `include/Rivet/Makefile.am` and `src/Projections/Makefile.am` files. These should be fairly easy to understand: just add your projection's header and implementation file names to the appropriate lists of projections to be built.

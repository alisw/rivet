## Physics tips and pitfalls

Rivet analyses can be very powerful if designed and constructed with the most general application in mind.
Choices that may seem arbitrary or irrelevant for your immediate application can enhance, or limit, the usefulness of your
routine to other users. Some issues to consider, based on experience, are discussed below.
There isn't always a "right" or "wrong" answer, but they are worth thinking about. In fact, many of these issues
are actually physics issues to do with making a measurement as unambiguous and model-independent, and hence durable, as possible.

### The final state

Think carefully about what is considered to be a "final state" particle for the purposes of input to your selection, including for example the starting point for jet finders. In general a final-state particle is something that could, in principle at least, interact with the detector. 

* Colour triplets (e.g. top, bottom quarks) are not final state particles. There are so many "top quark" measurements that Rivet will let you work around this, but be aware that calculations implemented in MC event generators do not guarantee the physicality of a "final state top". Obviously for lighter partons the situation is worse. Identifying heavy flavour hadrons is generally ok.

* Electroweak scale particles (W, Z, H) are not final state particles. Focus on the leptons and hadrons for the decay channel you care about.

* It may be useful to use lifetime cuts to define the final state, or to define a "prompt" final state in terms of pre-or-post hadronisation, or to use the rivet smearing utilities to get detector-level particles. All these things can be physically defensible and advantageous, depending on the analysis. 

[This document](https://cds.cern.ch/record/2022743?ln=en) contains some discussion which might be helpful.

### Hidden vetoes

Are all the important analysis cuts reflected in the fiducial cross section definition? For example, if a dilepton measurement vetoes on isolated photons, this may make no difference to the result when Rivet is run on a SM sample which is LO in electroweak corrections, but what happens if more precise calculation is used which may include EW radiation? And what happens if someone wants to evaluate possoble BSM contributions to the cross section? Ideally all important cuts should be reflected in the fiducial cross section definition, and thus implemented in your Rivet routine. If they are not in the fiducial measurement definition, you may want to nevertheless implement them approximately, and should at least warn users. 

### Missing energy and neutrinos

Many analyses "measure" neutrinos using missing transverse momentum, and may unfold to the true neutrino kinematics using a SM sample. This is not in itself a problem, but the Rivet routine should make use of the "true missing transverse momentum" projection rather than the true neutrino kinematics, so that possible additonal (BSM or other) sources of missing momentum can be correctly accounted for if present. This will give identical results if applied to the SM samples used for unfolding, but makes the routine more generally applicable.

Use of neutrino flavour and momentum explicitly in the measurement definition, especially in events with multiple neutrinos, is very problematic. While it is physically defined (modulo neutrino oscillatons!) it is of course not observable in practice and a measurement defined this way is certainly assuming the SM (there's no way of predicting how a DM candidate might contribute, for example). When trying to measure some very interesting final states (WZ processes etc) with as few assumptions as possible, a version of such measurements that can be implement in Rivet in terms of charged leptons and missing momentum is optimal. 





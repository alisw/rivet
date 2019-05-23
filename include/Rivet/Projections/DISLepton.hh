// -*- C++ -*-
#ifndef RIVET_DISLepton_HH
#define RIVET_DISLepton_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/UndressBeamLeptons.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Get the incoming and outgoing leptons in a DIS event.
  class DISLepton : public Projection {
  public:

    /// Enum to enable different orderings for selecting scattered
    /// leptons in case several were found.
    enum SortOrder { ENERGY, ETA, ET };
    
    /// @name Constructors.
    //@{

    /// Default constructor taking general options. The recognised
    /// options are: LMODE, taking the options "prompt", "any" and
    /// "dressed"; DressedDR giving a delta-R cone radius where photon
    /// momenta are added to the lepton candidates for LMODE=dresses;
    /// IsolDR giving a cone in delta-R where no hadrons are allowed
    /// around a lepton candidate; and Undress giving a cone around
    /// the incoming incoming beam in which photons are considered
    /// initial state rafiation for which the momentum is subtracted
    /// from the beam momentum.
    DISLepton(const std::map<std::string,std::string> & opts =
              std::map<std::string,std::string>())
      : _isolDR(0.0), _sort(ENERGY) {
      setName("DISLepton");
      addProjection(HadronicFinalState(), "IFS");

      auto sorting = opts.find("LSort");
      if ( sorting != opts.end() && sorting->second == "ETA" )
        _sort = ETA;
      else if ( sorting != opts.end() && sorting->second == "ET" )
        _sort = ET;

      double undresstheta = 0.0;
      auto undress = opts.find("Undress");
      if ( undress != opts.end() )
        undresstheta = std::stod(undress->second);
      if ( undresstheta > 0.0 )
        addProjection(UndressBeamLeptons(undresstheta), "Beam");
      else
        addProjection(Beam(), "Beam");

      auto isol = opts.find("IsolDR");
      if ( isol != opts.end() ) _isolDR = std::stod(isol->second);

      double dressdr = 0.0;
      auto dress = opts.find("DressDR");
      if ( dress != opts.end() )
        dressdr = std::stod(dress->second);

      auto lmode = opts.find("LMode");
      if ( lmode != opts.end() && lmode->second == "any" )
        addProjection(FinalState(), "LFS");
      else if ( lmode != opts.end() && lmode->second == "dressed" )
        addProjection(DressedLeptons(dressdr), "LFS");
      else
        addProjection(PromptFinalState(), "LFS");
    }

    /// Constructor taking the following arguments: a final state
    /// projection defining which lepton candidates to consider; a
    /// beam projection detining the momenta of the incoming lepton
    /// beam, and a final state projection defining the particles not
    /// allowed witin a delta-R of @a isolationcut of a lepton
    /// candidate.
    DISLepton(const FinalState & leptoncandidates,
              const Beam &  beamproj = Beam(),
              const FinalState & isolationfs = FinalState(),
              double isolationcut = 0.0, SortOrder sorting = ENERGY)
      : _isolDR(isolationcut), _sort(sorting) {
      addProjection(leptoncandidates, "LFS");
      addProjection(isolationfs, "IFS");
      addProjection(beamproj, "Beam");
    }


    
    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(DISLepton);

    //@}


  protected:

    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection& p) const;


  public:

    /// The incoming lepton
    const Particle& in() const { return _incoming; }

    /// The outgoing lepton
    const Particle& out() const { return _outgoing; }

    /// Sign of the incoming lepton pz component
    int pzSign() const { return sign(_incoming.pz()); }


  private:

    /// The incoming lepton
    Particle _incoming;

    /// The outgoing lepton
    Particle _outgoing;

    /// If larger than zerp an isolation cut around the lepton is required.
    double _isolDR;

    /// How to sort leptons
    SortOrder _sort;

  };

}


#endif

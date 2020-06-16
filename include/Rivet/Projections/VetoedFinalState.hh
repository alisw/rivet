// -*- C++ -*-
#ifndef RIVET_VetoedFinalState_HH
#define RIVET_VetoedFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief FS modifier to exclude classes of particles from the final state.
  class VetoedFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Constructor with a specific FinalState and a cuts list to veto
    VetoedFinalState(const FinalState& fsp, const vector<Cut>& cuts)
      : FinalState(), _vetoCuts(cuts)
    {
      setName("VetoedFinalState");
      declare(fsp, "FS");
    }

    /// Constructor with a specific FinalState and a single cut to veto
    VetoedFinalState(const FinalState& fsp, const Cut& cut)
      : VetoedFinalState(fsp, vector<Cut>{cut})
    {   }

    /// Constructor with a default FinalState and a cuts list to veto
    VetoedFinalState(const vector<Cut>& cuts)
      : VetoedFinalState(FinalState(), cuts)
    {   }

    /// Constructor with a default FinalState and a single cut to veto
    VetoedFinalState(const Cut& cut)
      : VetoedFinalState(FinalState(), vector<Cut>{cut})
    {   }

    /// Constructor with a specific FinalState and a PID list to veto
    VetoedFinalState(const FinalState& fsp, const vector<PdgId>& vetopids)
      : VetoedFinalState(fsp, {})
    {
      _vetoCuts.reserve(vetopids.size());
      for (PdgId pid : vetopids) addVeto(pid);
    }

    /// Constructor with a specific FinalState and a PID to veto
    VetoedFinalState(const FinalState& fsp, PdgId vetopid)
      : VetoedFinalState(fsp, vector<Cut>{Cuts::pid == vetopid})
    {   }

    /// Constructor with a default FinalState and a PID list to veto
    VetoedFinalState(const vector<PdgId>& vetopids)
      : VetoedFinalState(FinalState(), {})
    {
      _vetoCuts.reserve(vetopids.size());
      for (PdgId pid : vetopids) addVeto(pid);
    }

    /// Constructor with a default FinalState and a PID to veto
    VetoedFinalState(PdgId vetopid)
      : VetoedFinalState(FinalState(), vector<Cut>{Cuts::pid == vetopid})
    {   }

    /// Constructor with specific FinalState but no cuts
    VetoedFinalState(const FinalState& fsp)
      : VetoedFinalState(fsp, vector<Cut>())
    {   }

    /// Default constructor with default FinalState and no cuts
    VetoedFinalState()
      : VetoedFinalState(FinalState(), vector<Cut>())
    {   }

    /// You can add a map of ID plus a pair containing \f$ p_{Tmin} \f$ and
    /// \f$ p_{Tmax} \f$ -- these define the range of particles to be vetoed.
    //DEPRECATED("Prefer constructors using Cut arguments")
    VetoedFinalState(const map<PdgId,pair<double,double>>& vetocodes)
      : VetoedFinalState(FinalState(), {})
    {
      for (const auto& it : vetocodes) {
        addVeto(it.first, Cuts::pT > it.second.first && Cuts::pT < it.second.second);
      }
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(VetoedFinalState);

    //@}



    /// Get the list of particle IDs and \f$ p_T \f$ ranges to veto.
    const vector<Cut>& vetoDetails() const {
      return _vetoCuts;
    }
    //using vetos = vetoDetails;


    /// Add a particle selection to be vetoed from the final state
    VetoedFinalState& addVeto(const Cut& cut) {
      _vetoCuts.push_back(cut);
      return *this;
    }

    /// Add a particle selection to be vetoed from the final state
    VetoedFinalState& addVeto(PdgId pid, const Cut& cut=Cuts::OPEN) {
      _vetoCuts.push_back(Cuts::pid == pid && cut);
      return *this;
    }

    /// Add a particle/antiparticle selection to be vetoed from the final state
    VetoedFinalState& addVetoPair(PdgId pid, const Cut& cut=Cuts::OPEN) {
      _vetoCuts.push_back(Cuts::abspid == pid && cut);
      return *this;
    }


    /// @brief Add a particle ID and \f$ p_T \f$ range to veto
    ///
    /// Particles with \f$ p_T \f$ IN the given range will be rejected.
    VetoedFinalState& addVetoDetail(PdgId pid, double ptmin, double ptmax=std::numeric_limits<double>::max()) {
      return addVeto(pid, Cuts::ptIn(ptmin, ptmax));
    }
    //const auto addVeto = addVetoDetail;

    /// @brief Add a particle/antiparticle pair to veto in a given \f$ p_T \f$ range
    ///
    /// Given a single ID, both the particle and its conjugate antiparticle will
    /// be rejected if their \f$ p_T \f$ is IN the given range.
    VetoedFinalState& addVetoPairDetail(PdgId pid, double ptmin, double ptmax=std::numeric_limits<double>::max()) {
      return addVetoPair(pid, Cuts::ptIn(ptmin, ptmax));
    }
    //using addVetoPair = addVetoPairDetail;

    /// @brief Add a particle ID to veto (all \f$ p_T \f$ range will be vetoed)
    VetoedFinalState& addVetoId(PdgId pid) {
      return addVeto(pid);
    }
    //using addVeto = addVetoId;

    ///  @brief Add a particle/antiparticle pair to veto
    ///
    /// Given a single ID, both the particle and its corresponding antiparticle
    /// (for all \f$ p_T \f$ values) will be vetoed.
    VetoedFinalState& addVetoPairId(PdgId pid) {
      return addVetoPair(pid);
    }
    //using addVetoPair = addVetoPairId;


    /// Set the list of particle selections to veto
    VetoedFinalState& setVetoDetails(const vector<Cut>& cuts) {
      _vetoCuts = cuts;
      return *this;
    }
    //const auto setVetos = setVetoDetails;


    /// Veto all neutrinos (convenience method)
    VetoedFinalState& vetoNeutrinos() {
      addVetoPairId(PID::NU_E);
      addVetoPairId(PID::NU_MU);
      addVetoPairId(PID::NU_TAU);
      return *this;
    }


    /// Add a veto on composite masses within a given width.
    /// The composite mass is composed of nProducts decay products
    /// @ todo might we want to specify a range of pdg ids for the decay products?
    VetoedFinalState& addCompositeMassVeto(double mass, double width, int nProducts=2) {
      const double halfWidth = 0.5*width;
      pair<double,double> massRange(mass-halfWidth, mass+halfWidth);
      _compositeVetoes.insert(make_pair(nProducts, massRange));
      _nCompositeDecays.insert(nProducts);
      return *this;
    }


    /// Veto the decay products of particle with PDG code @a pid
    /// @todo Need HepMC to sort themselves out and keep vector bosons from
    /// the hard vtx in the event record before this will work reliably for all pdg ids
    VetoedFinalState& addDecayProductsVeto(PdgId pid) {
      _parentVetoes.insert(pid);
      return *this;
    }

    /// Veto particles from a supplied final state
    VetoedFinalState& addVetoOnThisFinalState(const ParticleFinder& fs) {
      const string name = "FS_" + to_str(_vetofsnames.size());
      declare(fs, name);
      _vetofsnames.insert(name);
      return *this;
    }


    /// Clear the list of particle IDs and ranges to veto.
    VetoedFinalState& reset() {
      _vetoCuts.clear();
      return *this;
    }


    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  private:

    /// The veto cuts
    vector<Cut> _vetoCuts;

    /// Composite particle masses to veto
    /// @todo Also generalise to use Cuts
    multimap<PdgId, pair<double,double> > _compositeVetoes;
    set<int> _nCompositeDecays;

    /// Set of decaying particle IDs to veto
    set<PdgId> _parentVetoes;

    /// Set of finalstate to be vetoed
    set<string> _vetofsnames;

  };


}


#endif

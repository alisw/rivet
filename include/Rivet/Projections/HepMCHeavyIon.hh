// -*- C++ -*-
#ifndef RIVET_HepMCHeavyIon_HH
#define RIVET_HepMCHeavyIon_HH

#include "Rivet/Projection.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  class HepMCHeavyIon : public Projection {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    HepMCHeavyIon();

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(HepMCHeavyIon);

    //@}


  protected:

    /// Perform the projection on the Event
    void project(const Event& e);

    /// Compare with other projections
    //int compare(const Projection& p) const;
    // Taken from Thrust.hh
    CmpState compare(const Projection& p) const {
      return CmpState::EQ;
    }

  public:

    /// Check that there were at all any heavy ion info in HepMC
    bool ok() const { return _hi != nullptr; }

    /// @brief the number of hard nucleon-nucleon collisions.
    ///
    /// Model-dependent. Usually the number of nucleon-nucleon
    /// collisions containing a special signal process. A negative
    /// value means that the information is not available.
    int    Ncoll_hard() const;
    
    /// @brief the number of participating nucleons in the projectile.
    ///
    /// The number of nucleons in the projectile participating in an
    /// inelastic collision (see Ncoll). A negative value means that
    /// the information is not available.
    int    Npart_proj() const;
    
    /// @brief the number of participating nucleons in the target.
    ///
    /// The number of nucleons in the target participating in an
    /// inelastic collision (see Ncoll). A negative value means that
    /// the information is not available.
    int    Npart_targ() const;

    /// @brief the number of inelastic nucleon-nucleon collisions.
    ///
    /// Note that a one participating nucleon can be involved in many
    /// inelastic collisions, and that inelastic also includes
    /// diffractive excitation. A negative value means that the
    /// information is not available.
    /// 
    int    Ncoll() const;
    /// @brief Collisions with a diffractively excited target nucleon.
    ///
    /// The number of single diffractive nucleon-nucleon collisions
    /// where the target nucleon is excited. A negative value means
    /// that the information is not available.
    int    N_Nwounded_collisions() const;
    
    /// @brief Collisions with a diffractively excited projectile nucleon.
    ///
    /// The number of single diffractive nucleon-nucleon collisions
    /// where the projectile nucleon is excited. A negative value
    /// means that the information is not available.
    int    Nwounded_N_collisions() const;

    /// @brief Non-diffractive or doubly diffractive collisions.
    ///
    /// The number of nucleon-nucleon collisions where both projectile
    /// and target nucleons are wounded. A negative value means that
    /// the information is not available.
    int    Nwounded_Nwounded_collisions() const;

    /// @brief The impact parameter.
    ///
    /// The impact parameter given in units of femtometer. A negative
    /// value means that the information is not available.
    double impact_parameter() const;

    /// @brief The event plane angle.
    ///
    /// The angle wrt. the x-axix of the impact parameter vector
    /// (pointing frm the target to the projectile). A positive number
    /// between 0 and two pi. A negative value means that the
    /// information is not available.
    double event_plane_angle() const;
    /// @brief The assumed inelastic nucleon-nucleon cross section
    ///
    /// in units of millibarn. As used in a Glauber calculation to
    /// simulate the distribution in Ncoll. A negative value means
    /// that the information is not available.
    double sigma_inel_NN() const;

    /// @brief The centrality.
    ///
    /// The generated centrality in percentiles, where 0 is the
    /// maximally central and 100 is the minimally central. A negative
    /// value means that the information is not available.
    double centrality() const;

    /// @brief A user defined centrality estimator.
    ///
    /// This variable may contain anything a generator feels is
    /// reasonable for estimating centrality. The value should be
    /// non-negative, and a low value corresponds to a low
    /// centrality. A negative value indicatess that the information
    /// is not available.
    double user_cent_estimate() const;

    /// @brief The number of spectator neutrons in the projectile
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_proj_n() const;

    /// @brief The number of spectator neutrons in the target
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_targ_n() const;

    /// @brief The number of spectator protons in the projectile
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_proj_p() const;

    /// @brief The number of spectator protons in the target
    ///
    /// ie. those that thave not participated in any inelastic
    /// nucleon-nucleon collision. A negative value indicatess that
    /// the information is not available.
    int Nspec_targ_p() const;

    /// @brief Participant plane angles
    ///
    /// calculated to different orders. The key of the map specifies
    /// the order, and the value gives to the angle wrt. the
    /// event plane.
    map<int,double> participant_plane_angles() const;

    /// @brief Eccentricities
    ///
    /// Calculated to different orders. The key of the map specifies
    /// the order, and the value gives the corresponding eccentricity.
    map<int,double> eccentricities() const;

  private:

    /// A pointer to the actual heavy ion object
    ConstGenHeavyIonPtr _hi;

  };
}


#endif

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"

namespace Rivet {


/// @brief Scaled momenta of identified particles
class ZEUS_2011_I945935 : public Analysis {
public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ZEUS_2011_I945935);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

        const UnstableParticles& labcut = UnstableParticles();
        declare(labcut, "UFS");
        const DISKinematics& diskin = DISKinematics();
        declare(diskin, "Kinematics");
        const DISFinalState&  disfsbf = DISFinalState(labcut, DISFinalState::BoostFrame::BREIT, diskin);
        declare(disfsbf, "FSBF");

        for (size_t offset = 0; offset < 5; offset++) {
            book(_h_K0S[offset],  2, 1, offset+1);
            book(_h_LAMBDA[offset], 5, 1, offset+1);
        }

        book(_h_K0S[5], 3, 1, 1);
        book(_h_K0S[6], 3, 1, 2);
        book(_h_LAMBDA[5], 6, 1, 1);
        book(_h_LAMBDA[6], 6, 1, 2);
        book(_h_Q2_tmp, "_TMP/N", 7, 0, 7);
    }

    int getbinQ2v1(const DISKinematics& dk) {
        if (inRange(dk.Q2()/GeV2, 10.0, 40.0) && inRange(dk.x(), 0.001, 0.75) ) return 1;
        if (inRange(dk.Q2()/GeV2, 40.0, 160.0) && inRange(dk.x(), 0.001, 0.75) ) return 2;
        if (inRange(dk.Q2()/GeV2, 160.0, 640.0) && inRange(dk.x(), 0.001, 0.75)) return 3;
        if (inRange(dk.Q2()/GeV2, 640.0, 2560.0) && inRange(dk.x(), 0.001, 0.75)) return 4;
        if (inRange(dk.Q2()/GeV2, 2650.0, 10240.0) && inRange(dk.x(), 0.001, 0.75)) return 5;
        return -1;
    }
    int getbinQ2v2(const DISKinematics& dk) {
        if (inRange(dk.Q2()/GeV2, 10.0, 100.0) && inRange(dk.x(), 0.001, 0.75) ) return 1;
        if (inRange(dk.Q2()/GeV2, 100.0, 1000.0) && inRange(dk.x(), 0.001, 0.75) ) return 2;
        return -1;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
        /// DIS kinematics
        const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
        const double q2  = dk.Q2();
        const double x = dk.x();
        const double y = dk.y();

        if (!inRange(q2/GeV2, 10.0, 40000.0)) vetoEvent;
        if (!inRange(y, 0.04, 0.95)) vetoEvent;
        if (!inRange(x, 0.001, 0.75)) vetoEvent;
        const int ofv1 = getbinQ2v1(dk) - 1;
        const int ofv2 = getbinQ2v2(dk) - 1;
        if ( ofv1 < 0 ) vetoEvent;

        /// @todo Do the event by event analysis here
        _h_Q2_tmp->fill(ofv1);
        if (ofv2 >= 0) _h_Q2_tmp->fill(5+ofv2);
        const DISFinalState& disfsbf = apply<DISFinalState>(event, "FSBF");

        for (const Particle& p: filter_select(disfsbf.particles(), Cuts::abspid == abs(PID::K0S))) {
            //// Scaled energy.
            if (p.pz() > 0) continue;
            const double energy = p.momentum().vector3().mod();
            const double scaledEnergy = 2.0*energy/sqrt(q2);
            _h_K0S[ofv1]->fill(scaledEnergy);
            if (ofv2 >= 0) _h_K0S[5+ofv2]->fill(scaledEnergy);
        }

        for (const Particle& p: filter_select(disfsbf.particles(), Cuts::abspid == abs(PID::LAMBDA))) {
            //// Scaled energy.
            if (p.pz() > 0) continue;
            const double energy = p.momentum().vector3().mod();
            const double scaledEnergy = 2.0*energy/sqrt(q2);
            _h_LAMBDA[ofv1]->fill(scaledEnergy);
            if (ofv2 >= 0) _h_LAMBDA[5+ofv2]->fill(scaledEnergy);
        }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
        for (size_t offset = 0; offset < 7; offset++) {
            scale(_h_K0S[offset], 1.0/_h_Q2_tmp->bin(offset).area());
            scale(_h_LAMBDA[offset], 1.0/_h_Q2_tmp->bin(offset).area());
        }
    }
    Histo1DPtr _h_K0S[7];
    Histo1DPtr _h_LAMBDA[7];
    Histo1DPtr _h_Q2_tmp;
};
/// The hook for the plugin system
DECLARE_RIVET_PLUGIN(ZEUS_2011_I945935);
}

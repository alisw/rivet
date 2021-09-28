// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {
    
    /// @brief single-diffractive cross-sections at 8 TeV
    class ATLAS_2019_I1762584 : public Analysis {
        public:
        
        /// Constructor
        DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2019_I1762584);
        
        void init() {
            
            // Inner Detector (ID) charged final states
            // charged tracks with pT>200MeV, |eta|<2.5
            Cut track_cuts = Cuts::abseta < 2.5 && Cuts::pT > 0.2*GeV;
            const ChargedFinalState tracks(track_cuts);
            declare(tracks, "tracks");
            
            // Forward Detector (FD) protons
            // protons with 0.016<|t|/GeV^2<0.43: 0.126<pT/GeV<0.655
            Cut proton_t_cuts = Cuts::pT>0.126/GeV && Cuts::pT<0.655/GeV  ;
            // protons with -4<log10(Xi)<-1.6: 0.016<|t|/GeV^2<0.43: 3899.52<E<3999.6
            Cut proton_xi_cuts = Cuts::E>3899.52/GeV && Cuts::E<3999.6/GeV ;
            const ChargedFinalState protons(Cuts::pid==PID::PROTON && proton_t_cuts && proton_xi_cuts ) ;
            declare(protons, "protons");
            
            // Book histograms
            book(_h_dSigma_dDeltaEta, 1, 1, 1);
            book(_h_dSigma_dAbsT, 2, 1, 1);
            book(_h_dSigma_dLog10Xi, 3, 1, 1);
            
        }
        
        
        void analyze(const Event& event) {
            
            // Retrieve charged tracks in Inner Detector
            const ChargedFinalState& tracks = apply<ChargedFinalState>(event, "tracks");
            
            // Retrieve protons in ALFA detector
            const ChargedFinalState protons = apply<ChargedFinalState>(event, "protons");
            
            // Veto Events with more than one tagged proton
            if (protons.size()!=1) vetoEvent;
            const Particle tagProton = protons.particles()[0];
            
            // Calculate |t|
            const double AbsT = (tagProton.pT()/GeV)*(tagProton.pT()/GeV);
            
            // Calculate Log10(Xi)
            const double Log10Xi = log10(1.-(tagProton.E()/GeV)/4000.);
            
            //Calculate DeltaEta
            const double EtaEdge = 2.5*tagProton.pz()/abs(tagProton.pz());
            double DeltaEta = 5;
            for (const Particle& p : tracks.particles()) {
                double DeltaEta_track = abs(p.eta() - EtaEdge);
                if (DeltaEta_track<DeltaEta) DeltaEta=DeltaEta_track;
            }
            
            // Fill histograms
            _h_dSigma_dDeltaEta->fill(DeltaEta);
            _h_dSigma_dAbsT->fill(AbsT);
            _h_dSigma_dLog10Xi->fill(Log10Xi);
            
        }
        
        /// Normalise histograms to units of millibarn
        void finalize() {
            
            // norm to generated cross-section in mb (after cuts)
            scale(_h_dSigma_dAbsT, crossSection()/millibarn/sumOfWeights());
            scale(_h_dSigma_dLog10Xi, crossSection()/millibarn/sumOfWeights());
            scale(_h_dSigma_dDeltaEta, crossSection()/millibarn/sumOfWeights());
        }
        
        
    private:
        
        Histo1DPtr _h_dSigma_dAbsT, _h_dSigma_dLog10Xi, _h_dSigma_dDeltaEta;
        
    };
    
    DECLARE_RIVET_PLUGIN(ATLAS_2019_I1762584);
    
}


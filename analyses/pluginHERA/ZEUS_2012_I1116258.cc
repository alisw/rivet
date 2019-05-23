// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/SISConePlugin.hh"

namespace Rivet {


  /// @brief ZEUS inclusive jet photoproduction study used to measure alpha_s
  ///
  /// @author Jon Butterworth
  class ZEUS_2012_I1116258 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ZEUS_2012_I1116258);

    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {

      // Projections

      // Jet schemes checked with oringal code, M.Wing, A.Geiser
      FinalState fs;
      double jet_radius = 1.0;
      declare(FastJets(fs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::Et_scheme, jet_radius), "Jets"); 
      declare(FastJets(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::Et_scheme, jet_radius), "Jets_akt"); 

      // bit of messing about to use the correct recombnation scheme for SISCone.
      double overlap_threshold = 0.75;
      fastjet::SISConePlugin * plugin = new fastjet::SISConePlugin(jet_radius, overlap_threshold);
      plugin->set_use_jet_def_recombiner(true);
      JetDefinition siscone(plugin);
      siscone.set_recombination_scheme(fastjet::RecombinationScheme::Et_scheme);
      declare(FastJets(fs, siscone), "Jets_sis"); 

      
      declare(DISKinematics(), "Kinematics");
      
      // all eta
      _h_etjet[0] = bookHisto1D(1, 1, 1);
      
      // two ET cuts.
      _h_etajet[0] = bookHisto1D(2, 1, 1);
      _h_etajet[1] = bookHisto1D(3, 1, 1);
      
      // in eta regions
      _h_etjet[1] = bookHisto1D(4, 1, 1);
      _h_etjet[2] = bookHisto1D(5, 1, 1);
      _h_etjet[3] = bookHisto1D(6, 1, 1);
      _h_etjet[4] = bookHisto1D(7, 1, 1);
      _h_etjet[5] = bookHisto1D(8, 1, 1);
      
      // antiKT
      _h_etjet[6] = bookHisto1D(9, 1, 1);
      _h_etajet[2] = bookHisto1D(11, 1, 1);
      
      // SiSCone
      _h_etjet[7] = bookHisto1D(10, 1, 1);
      _h_etajet[3] = bookHisto1D(12, 1, 1);
      
    }
    
    
    // Do the analysis
    void analyze(const Event& event) {
      
      // Determine kinematics, including event orientation since ZEUS coord system is for +z = proton direction
      const DISKinematics& kin = apply<DISKinematics>(event, "Kinematics");
      const int orientation = kin.orientation();
      
      // Q2 and inelasticity cuts
      if (kin.Q2() > 1*GeV2) vetoEvent;
      if (!inRange(sqrt(kin.W2()), 142.0, 293.0)) vetoEvent;
      
      // Jet selection
      // @TODO check the recombination scheme
      const Jets jets = apply<FastJets>(event, "Jets")			\
        .jets(Cuts::Et > 17*GeV && Cuts::etaIn(-1*orientation, 2.5*orientation), cmpMomByEt);
      MSG_DEBUG("kT Jet multiplicity = " << jets.size());
      
      const Jets jets_akt = apply<FastJets>(event, "Jets_akt")		\
        .jets(Cuts::Et > 17*GeV && Cuts::etaIn(-1*orientation, 2.5*orientation), cmpMomByEt);

      const Jets jets_sis = apply<FastJets>(event, "Jets_sis")	\
	.jets(Cuts::Et > 17*GeV && Cuts::etaIn(-1*orientation, 2.5*orientation), cmpMomByEt);


      // Fill histograms
      const double weight = event.weight();

      for (const Jet& jet : jets ){
	_h_etjet[0]->fill(jet.pt(), weight);
	_h_etajet[0]->fill(orientation*jet.eta(), weight);
	if (jet.pt()>21*GeV) {
	  _h_etajet[1]->fill(orientation*jet.eta(), weight);
	}
	if (orientation*jet.eta() < 0) { 
	  _h_etjet[1]->fill(jet.pt(), weight);
	} else if (orientation*jet.eta() < 1) { 
	  _h_etjet[2]->fill(jet.pt(), weight);
	} else if (orientation*jet.eta() < 1.5) { 
	  _h_etjet[3]->fill(jet.pt(), weight);
	} else if (orientation*jet.eta() < 2) { 
	  _h_etjet[4]->fill(jet.pt(), weight);
	} else { 
	  _h_etjet[5]->fill(jet.pt(), weight);
	}
      }

      for (const Jet& jet : jets_akt ){
	_h_etjet[6]->fill(jet.pt(), weight);
	_h_etajet[2]->fill(orientation*jet.eta(), weight);
      }
      for (const Jet& jet : jets_sis ){
	_h_etjet[7]->fill(jet.pt(), weight);
	_h_etajet[3]->fill(orientation*jet.eta(), weight);
      }

    }


      // Finalize
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();
      for( int i = 0; i < 8; i++ ) {
	scale(_h_etjet[i], sf);
      }
      for( int i = 0; i < 4; i++ ) {
	scale(_h_etajet[i], sf);
      }
    }

    //@}
    

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_etjet[8], _h_etajet[4]; 
    //@}

    };
  
  
  DECLARE_RIVET_PLUGIN(ZEUS_2012_I1116258);
  
}

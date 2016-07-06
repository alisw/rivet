/*
 * Rivet routine for the Higgs differential cross sections combination
 * between the ATLAS measurements in the yy and 4l channels for
 * Higgs transverse momentum, rapidity, jet multiplicity and leading jet pT
 *
 * Author: Michaela Queitsch-Maitland (ATLAS Collaboration)
 * Contact: Michaela Queitsch-Maitland <michaela.queitsch-maitland@cern.ch>,
 *          Dag Gillberg <dag.gillberg@cern.ch>,
 *          Florian Bernlochner <florian.bernlochner@cern.ch>,
 *          Sarah Heim <sarah.heim@cern.ch>
 */

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vector4.hh"
#include "Rivet/Particle.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {


  class ATLAS_2015_I1364361 : public Analysis {
  public:

    /// Constructor
    ATLAS_2015_I1364361()
      : Analysis("ATLAS_2015_I1364361")
    {    }



    /// Book histograms and initialise projections before the run
    void init() {

      // All final state particles
      const FinalState FS;
      addProjection(FS,"FS");

      // Histograms with data bins
      _h_pTH_incl   = bookHisto1D(1,1,1);  
      _h_yH_incl    = bookHisto1D(2,1,1);   
      _h_Njets_incl = bookHisto1D(3,1,1);
      _h_pTj1_incl  = bookHisto1D(4,1,1); 
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get event weight
      const double weight = event.weight();

      // Get the final state particles ordered by pT
      const ParticleVector& FS = applyProjection<FinalState>(event, "FS").particlesByPt();

      // Find the Higgs
      bool stable_higgs = false;
      const Particle* higgs=0;
      foreach ( const Particle& p, FS ) {
	if ( p.pid()==25 ) {
	  stable_higgs = true;
	  higgs = &p;
	  break;
	}
      }

      // If no stable Higgs found in event record, can't do anything
      if ( !stable_higgs ) {
	MSG_WARNING("FATAL: No stable Higgs found in event record.\n");
	vetoEvent;
      }

      ParticleVector leptons;
      ParticleVector photons; // for dressing
      ParticleVector jet_ptcls;

      // Loop over final state particles and fill jet particles vector
      foreach ( const Particle& ptcl, FS ) {
	// Do not include the Higgs in jet finding!
	if ( ptcl.pid()==25 ) continue;
	// Neutrinos not from hadronisation
	if ( ptcl.isNeutrino() && !ptcl.fromHadron() ) continue;
	// Electrons and muons not from hadronisation
	if ( ( ptcl.abspid() == 11 || ptcl.abspid() == 13 ) && !ptcl.fromHadron() ) {
	  leptons.push_back(ptcl);
	  continue;
	}
	// Photons not from hadronisation
	if ( ptcl.abspid() == 22 && !ptcl.fromHadron() ) {
	  photons.push_back(ptcl);
	  continue;
	}
	// Add particle to jet inputs
	jet_ptcls.push_back(ptcl);
      }

      // Match FS photons to leptons within cone R=0.1
      // If they are not 'dressing' photons, add to jet particle vector
      foreach ( const Particle& ph, photons ) {
	bool fsr_photon = false;
	foreach ( const Particle& lep, leptons ) {
	  if ( deltaR(ph.momentum(),lep.momentum()) < 0.1 ){
	    fsr_photon=true;
	    continue;
	  }
	}
	if ( !fsr_photon ) jet_ptcls.push_back(ph);
      }

      // Let's build the jets!
      FastJets jet_pro(FastJets::ANTIKT, 0.4);
      jet_pro.calc(jet_ptcls);
      Jets jets = jet_pro.jetsByPt(Cuts::pT>30*GeV && Cuts::absrap<4.4);

      _pTH = higgs->momentum().pT();
      _yH = higgs->momentum().rapidity();
      _Njets = jets.size() > 3 ? 3 : jets.size();
      _pTj1 = jets.size() > 0 ? jets[0].momentum().pT() : 0.;

      _h_pTH_incl->fill(_pTH,weight);
      _h_yH_incl->fill( fabs(_yH),weight);
      _h_Njets_incl->fill(_Njets,weight);
      _h_pTj1_incl->fill(_pTj1,weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double xs = crossSectionPerEvent();
      scale( _h_pTH_incl, xs );
      scale( _h_yH_incl, xs );
      scale( _h_Njets_incl, xs );
      scale( _h_pTj1_incl, xs );
    }



  private:

    double _pTH;
    double _yH;
    double _Njets;
    double _pTj1;


    Histo1DPtr _h_pTH_incl;
    Histo1DPtr _h_yH_incl;
    Histo1DPtr _h_Njets_incl;
    Histo1DPtr _h_pTj1_incl;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1364361);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

/// @brief Measurement of forward neutron(+ antineutron) production in proton-proton Collisions at 13 TeV
class LHCF_2020_I1783943 : public Analysis {
public:

	/// Constructor
	RIVET_DEFAULT_ANALYSIS_CTOR(LHCF_2020_I1783943);

	/// @name Analysis methods
	//@{

	/// Book histograms and initialise projections before the run
	void init() {

		// Initialise and register projections
		declare(Beam(), "Beam");
		declare(FinalState(), "FS");

		// Book histograms
		book(_h_n_en_eta1, 1, 1, 1);
		book(_h_n_en_eta2, 2, 1, 1);
		book(_h_n_en_eta3, 3, 1, 1);
		book(_h_n_en_eta4, 4, 1, 1);
		book(_h_n_en_eta5, 5, 1, 1);
		book(_h_n_en_eta6, 6, 1, 1);
		book(_h_n_eflow, 7, 1, 1);
		book(_h_n_sigma, 8, 1, 1);
		book(_h_n_elas,  9, 1, 1);
		book(_h_n_inel, 10, 1, 1);

		book(_inelnorm, "_inelasticity_norm");
	}

	/// Perform the per-event analysis
	void analyze(const Event& event) {

		const FinalState &fs = apply<FinalState> (event, "FS");
		Particles fs_particles = fs.particles();

		bool is_fl_neutron = false;
		double elasticity = 0.;

		for (Particle& p: fs_particles ) {
			//Artificially remove QGSJet II-04 wrong events
			//if(p.pT()/GeV == 0.0) continue;

			// double analysis efficiency with a two-sided LHCf
			double energy = p.E()/GeV;
			double eta = abs(p.eta());
			// highest eta value must be in [10.75, 13.00]
			if(eta > 13.0)
				eta = 11.875;

			// search for forward leading particle on one side only
			if(p.eta() > 0) {
				if(2.*energy/sqrtS() > elasticity) {
					elasticity = 2.*energy/sqrtS();
					if(p.abspid() == 2112)
						is_fl_neutron = true;
					else
						is_fl_neutron = false;
				}
			}

			// select neutrons above threshold
			if(p.abspid() != 2112) continue; //Select neutrons and antineutrons only

			_h_n_eflow->fill( eta , energy ); //Energy Flow
			_h_n_sigma->fill( eta ,    1.0 ); //Cross Section

			// select only energy above 500 GeV
			if(p.E()/GeV < 500.) continue;

			// fill energy distributions
			if( eta >= 10.749356 ) {
				_h_n_en_eta1->fill( energy );
			} else if(eta >= 10.056209 && eta <= 10.749356) {
				_h_n_en_eta2->fill( energy );
			} else if(eta >= 9.650744 && eta <= 10.056209) {
				_h_n_en_eta3->fill( energy );
			} else if(eta >= 8.985767 && eta <= 9.208911) {
				_h_n_en_eta4->fill( energy );
			} else if(eta >= 8.803446 && eta <= 8.985767) {
				_h_n_en_eta5->fill( energy );
			} else if(eta >= 8.649295 && eta <= 8.803446) {
				_h_n_en_eta6->fill( energy );
			}

		}

		// fill elasticity distribution if forward particle is a neutron
		if(is_fl_neutron) {
			_h_n_elas->fill(         elasticity );
			_h_n_inel->fill( 1.0, 1.-elasticity );
			_inelnorm->fill();
		}

	}

	/// Normalise histograms etc., after the run
	void finalize() {

		//Scale considering the LHCf Arm2 side
		scale(_h_n_en_eta1, crossSection()/millibarn/sumOfWeights()/2.); // normalize to cross section
		scale(_h_n_en_eta2, crossSection()/millibarn/sumOfWeights()/2.); // normalize to cross section
		scale(_h_n_en_eta3, crossSection()/millibarn/sumOfWeights()/2.); // normalize to cross section
		scale(_h_n_en_eta4, crossSection()/millibarn/sumOfWeights()/2.); // normalize to cross section
		scale(_h_n_en_eta5, crossSection()/millibarn/sumOfWeights()/2.); // normalize to cross section
		scale(_h_n_en_eta6, crossSection()/millibarn/sumOfWeights()/2.); // normalize to cross section

		scale(_h_n_eflow, 1./sumOfWeights()/2.); // normalize to event number
		scale(_h_n_sigma, crossSection()/millibarn/sumOfWeights()/2.); // normalize to cross section

		scale(_h_n_inel, 1. / *_inelnorm); // normalize to forward leading neutron event number
		for (unsigned int ibin=0; ibin<_h_n_inel->numBins(); ++ibin) // multiply for bin width
			_h_n_inel->bin(ibin).scaleW(_h_n_inel->bin(ibin).xMax()-_h_n_inel->bin(ibin).xMin());
		scale(_h_n_elas, crossSection()/millibarn/sumOfWeights()); // normalize to cross section
		for (unsigned int ibin=0; ibin<_h_n_elas->numBins(); ++ibin) // multiply for bin width
			_h_n_elas->bin(ibin).scaleW(_h_n_elas->bin(ibin).xMax()-_h_n_elas->bin(ibin).xMin());

	}
	//@}

private:

	/// @name Histograms
	//@{
	Histo1DPtr _h_n_en_eta1;
	Histo1DPtr _h_n_en_eta2;
	Histo1DPtr _h_n_en_eta3;
	Histo1DPtr _h_n_en_eta4;
	Histo1DPtr _h_n_en_eta5;
	Histo1DPtr _h_n_en_eta6;
	Histo1DPtr _h_n_eflow;
	Histo1DPtr _h_n_sigma;
	Histo1DPtr _h_n_inel;
	Histo1DPtr _h_n_elas;

	CounterPtr _inelnorm;
	//@}
};

// The hook for the plugin system
RIVET_DECLARE_PLUGIN(LHCF_2020_I1783943);

}

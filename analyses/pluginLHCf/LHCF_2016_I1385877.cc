// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"
namespace Rivet {


/// @brief Add a short analysis description here
class LHCF_2016_I1385877 : public Analysis {
public:

	/// Constructor
	DEFAULT_RIVET_ANALYSIS_CTOR(LHCF_2016_I1385877);

	//In case of some models there can be very small value pT but greater than 0.
	//In order to avoid unphysical behavior in the first bin a cutoff is needed
	//If you are sure the model does not have this problem you can set pt_cutoff to 0.
	const double pt_cutoff = 0.01;

	/// @name Analysis methods
	//@{

	/// Book histograms and initialise projections before the run
	void init() {

		// Initialise and register projections
		addProjection(UnstableParticles(), "UFS");
		addProjection(Beam(), "Beam");

		// calculate beam rapidity
		const Particle bm1 = beams().first;
		const Particle bm2 = beams().second;
		_beam_rap_1 = bm1.rap();
		_beam_rap_2 = bm2.rap();
		MSG_INFO("Beam 1 : momentum " << bm1.pz() << " PID " << bm1.pid() << " rapidity " << bm1.rap() );
		MSG_INFO("Beam 2 : momentum " << bm2.pz() << " PID " << bm2.pid() << " rapidity " << bm2.rap() );
		const double _sqrts = sqrtS();
		MSG_INFO("CM energy: " << _sqrts );

		_beam_rap = _beam_rap_1;
		if(bm1.pid()==2212 && bm2.pid()==2212) { //p-p
			_pp_Pb = true;
			if( fuzzyEquals( _sqrts/GeV, 7000., 1E-3) ) {
				_p_pi0_rap_apT = bookProfile1D(1, 1, 2);

				_h_pi0_rap_pT.addHistogram(  8.8, 9.0, bookHisto1D(2, 1, 2));
				_h_pi0_rap_pT.addHistogram(  9.0, 9.2, bookHisto1D(3, 1, 2));
				_h_pi0_rap_pT.addHistogram(  9.2, 9.4, bookHisto1D(4, 1, 2));
				_h_pi0_rap_pT.addHistogram(  9.4, 9.6, bookHisto1D(5, 1, 2));
				_h_pi0_rap_pT.addHistogram(  9.6, 9.8, bookHisto1D(6, 1, 2));
				_h_pi0_rap_pT.addHistogram(  9.8, 10.0, bookHisto1D(7, 1, 2));
				_h_pi0_rap_pT.addHistogram(  10.0, 10.2, bookHisto1D(8, 1, 2));
				_h_pi0_rap_pT.addHistogram(  10.2, 10.4, bookHisto1D(9, 1, 2));
				_h_pi0_rap_pT.addHistogram(  10.4, 10.6, bookHisto1D(10, 1, 2));
				_h_pi0_rap_pT.addHistogram(  10.6, 10.8, bookHisto1D(11, 1, 2));

				_h_pi0_pT_pZ.addHistogram(  0.0, 0.2, bookHisto1D(12, 1, 2));
				_h_pi0_pT_pZ.addHistogram(  0.2, 0.4, bookHisto1D(13, 1, 2));
				_h_pi0_pT_pZ.addHistogram(  0.4, 0.6, bookHisto1D(14, 1, 2));
				_h_pi0_pT_pZ.addHistogram(  0.6, 0.8, bookHisto1D(15, 1, 2));
				_h_pi0_pT_pZ.addHistogram(  0.8, 1.0, bookHisto1D(16, 1, 2));

				_h_pi0_rap = bookHisto1D(21, 1, 2);

				_p_pi0_raploss_apT = bookProfile1D(22, 1, 2);
				_h_pi0_raploss = bookHisto1D(23, 1, 2);
			} else if(fuzzyEquals( _sqrts/GeV, 2760., 1E-3) ){
				_p_pi0_rap_apT = bookProfile1D(1, 1, 1);

				_h_pi0_rap_pT.addHistogram(  8.8, 9.0, bookHisto1D(2, 1, 1));
				_h_pi0_rap_pT.addHistogram(  9.0, 9.2, bookHisto1D(3, 1, 1));
				_h_pi0_rap_pT.addHistogram(  9.2, 9.4, bookHisto1D(4, 1, 1));
				_h_pi0_rap_pT.addHistogram(  9.4, 9.6, bookHisto1D(5, 1, 1));
				_h_pi0_rap_pT.addHistogram(  9.6, 9.8, bookHisto1D(6, 1, 1));

				_h_pi0_pT_pZ.addHistogram(  0.0, 0.2, bookHisto1D(12, 1, 1));
				_h_pi0_pT_pZ.addHistogram(  0.2, 0.4, bookHisto1D(13, 1, 1));

				_h_pi0_rap = bookHisto1D(21, 1, 1);

				_p_pi0_raploss_apT = bookProfile1D(22, 1, 1);
				_h_pi0_raploss = bookHisto1D(23, 1, 1);
			}else{
				MSG_INFO("p-p collisions : energy out of range!");
			}
		} else if (bm1.pid()==2212 && bm2.pid()==1000822080){ //p-Pb
			_pp_Pb = false;
			if( fuzzyEquals( _sqrts/sqrt(208.)/GeV, 5020., 1E-3) ) {
				_p_pi0_rap_apT = bookProfile1D(1, 1, 3);

				_h_pi0_rap_pT.addHistogram(  8.8, 9.0, bookHisto1D(2, 1, 3));
				_h_pi0_rap_pT.addHistogram(  9.0, 9.2, bookHisto1D(3, 1, 3));
				_h_pi0_rap_pT.addHistogram(  9.2, 9.4, bookHisto1D(4, 1, 3));
				_h_pi0_rap_pT.addHistogram(  9.4, 9.6, bookHisto1D(5, 1, 3));
				_h_pi0_rap_pT.addHistogram(  9.6, 9.8, bookHisto1D(6, 1, 3));
				_h_pi0_rap_pT.addHistogram(  9.8, 10.0, bookHisto1D(7, 1, 3));
				_h_pi0_rap_pT.addHistogram(  10.0, 10.2, bookHisto1D(8, 1, 3));
				_h_pi0_rap_pT.addHistogram(  10.2, 10.4, bookHisto1D(9, 1, 3));
				_h_pi0_rap_pT.addHistogram(  10.4, 10.6, bookHisto1D(10, 1, 3));
				_h_pi0_rap_pT.addHistogram(  10.6, 10.8, bookHisto1D(11, 1, 3));

				_h_pi0_pT_pZ.addHistogram(  0.0, 0.2, bookHisto1D(12, 1, 3));
				_h_pi0_pT_pZ.addHistogram(  0.2, 0.4, bookHisto1D(13, 1, 3));
				_h_pi0_pT_pZ.addHistogram(  0.4, 0.6, bookHisto1D(14, 1, 3));
				_h_pi0_pT_pZ.addHistogram(  0.6, 0.8, bookHisto1D(15, 1, 3));
				_h_pi0_pT_pZ.addHistogram(  0.8, 1.0, bookHisto1D(16, 1, 3));

				//_h_pi0_rap = bookHisto1D(21, 1, 3);

				_p_pi0_raploss_apT = bookProfile1D(22, 1, 3);
				//_h_pi0_raploss = bookHisto1D(23, 1, 3);
			}else{
				MSG_INFO("p-Pb collisions : energy out of range!");
			}
		} else {
			MSG_INFO("Beam PDGID out of range!");
		}
		_nevt = 0.;
	}

	/// Perform the per-event analysis
	void analyze(const Event& event) {
		_nevt = _nevt + 1.;

		const UnstableParticles &ufs = applyProjection<UnstableFinalState> (event, "UFS");
		Particles ufs_particles = ufs.particles();
		for (Particle& p: ufs_particles ) {

			// select neutral pion
			if(p.abspid() != 111) continue;
			if(p.pz()/GeV<0.) continue;
			if(p.pT()/GeV<pt_cutoff) continue;

			// Kinematic variables
			if(_pp_Pb) {//p-p collisions
				const double pZ = p.pz()/GeV;
				const double pT = p.pT()/GeV;
				const double pT_MeV = p.pT()/MeV;
				const double en = p.E()/GeV;
				const double rap = p.rap();
				const double raploss = _beam_rap_1 - p.rap();

				_p_pi0_rap_apT->fill( rap , pT_MeV , 1.0 );
				_h_pi0_rap_pT.fill( rap, pT , 1.0 / pT );
				_h_pi0_pT_pZ.fill( pT, pZ , en / pT);
				_h_pi0_rap->fill( rap, 1.0 );
				_p_pi0_raploss_apT->fill( raploss , pT_MeV , 1.0 );
				_h_pi0_raploss->fill( raploss, 1.0 );
			} else {//pPb collisions
				const double pZ = p.pz()/GeV;
				const double pT = p.pT()/GeV;
				const double pT_MeV = p.pT()/MeV;
				const double en = p.E()/GeV;
				const double rap = p.rap();
				const double raploss = _beam_rap_1 - p.rap();  //mitsuka-like

				_p_pi0_rap_apT->fill( rap , pT_MeV , 1.0 );
				_h_pi0_rap_pT.fill( rap, pT , 1.0 / pT );
				_h_pi0_pT_pZ.fill( pT, pZ , en / pT);
				//_h_pi0_rap->fill( rap, 1.0 );
				_p_pi0_raploss_apT->fill( raploss , pT_MeV , 1.0 );
				//_h_pi0_raploss->fill( raploss, 1.0 );
			}

		}

	}

	/// Normalise histograms etc., after the run
	void finalize() {
		const double inv_scale_factor = 1. / _nevt / (2.*PI);
		const double pt_bin_width = 0.2;

		for (Histo1DPtr h: _h_pi0_pT_pZ.getHistograms()){
			if(h->path() == "/LHCF_2016_I1385877/d12-x01-y01" ||
					h->path() == "/LHCF_2016_I1385877/d12-x01-y02" ||
					h->path() == "/LHCF_2016_I1385877/d12-x01-y03")
				h->scaleW( inv_scale_factor / (pt_bin_width-pt_cutoff) );
			else
				h->scaleW( inv_scale_factor / pt_bin_width );
		}

		const double scale_factor =  1. / _nevt / (2.*PI);
		const double rap_bin_width = 0.2;
		for (Histo1DPtr h: _h_pi0_rap_pT.getHistograms()) {
			const int cutoff_bin = h->binIndexAt(pt_cutoff);

			if(cutoff_bin>=0) {
//				for(unsigned int ibin=0; ibin<h->numBins(); ++ibin)
//						cout << ibin << " " << h->bin(ibin).area() << endl;
				const double cutoff_wdt = h->bin(cutoff_bin).xMax()-h->bin(cutoff_bin).xMin();
				h->bin(cutoff_bin).scaleW((cutoff_wdt)/(cutoff_wdt-pt_cutoff));
//				for(unsigned int ibin=0; ibin<h->numBins(); ++ibin)
//						cout << ibin << " " << h->bin(ibin).area() << endl;
			}

			h->scaleW( scale_factor / rap_bin_width );
		}

		if(_pp_Pb) {
			scale( _h_pi0_rap , 1. / _nevt );
			scale( _h_pi0_raploss , 1. / _nevt );
		}
	}
	//@}


private:

	/// @name Histograms
	//@{
	bool _pp_Pb;
	double _nevt;
	double _beam_rap;
	double _beam_rap_1;
	double _beam_rap_2;
	BinnedHistogram<double> _h_pi0_pT_pZ;
	BinnedHistogram<double> _h_pi0_rap_pT;
	Profile1DPtr _p_pi0_rap_apT;
	Histo1DPtr _h_pi0_rap;
	Profile1DPtr _p_pi0_raploss_apT;
	Histo1DPtr _h_pi0_raploss;
	//@}
};

// The hook for the plugin system
DECLARE_RIVET_PLUGIN(LHCF_2016_I1385877);

}

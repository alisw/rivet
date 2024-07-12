// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
	// @brief ZEUS dijet photoproduction in different xgamma regions
  
	class ZEUS_1997_I450085 : public Analysis {
	public:

		// Constructor
		RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1997_I450085);
		
		// Initialization, books projections and histograms
		void init() {
			
			// Projections
			// checking recombination scheme and radius checked with original code from M.Wing
			FinalState fs;
			declare(FastJets(fs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::Et_scheme, 1.0), "Jets");
			declare(DISKinematics(), "Kinematics");

			//Histograms
			// Table 1
			book(_h_etabar_all[0], 1, 1, 1);
			book(_h_etabar_all[1], 2, 1, 1);
			book(_h_etabar_all[2], 3, 1, 1);
			book(_h_etabar_all[3], 4, 1, 1);
			// Table 6
			book(_h_etabar[1][0], 21, 1, 1);
			book(_h_etabar[1][1], 22, 1, 1);
			book(_h_etabar[1][2], 23, 1, 1);
			book(_h_etabar[1][3], 24, 1, 1);
			// Table 7
			book(_h_etabar[0][0], 25, 1, 1);
			book(_h_etabar[0][1], 26, 1, 1);
			book(_h_etabar[0][2], 27, 1, 1);
			book(_h_etabar[0][3], 28, 1, 1);

		}
		// Analysis
		void analyze(const Event & event) {

			// Determine kinematics, including event orientation since ZEUS coord system is for +z = proton direction
			const DISKinematics & kin = apply<DISKinematics>(event, "Kinematics");
		    if (kin.failed()) vetoEvent;
			const int orientation = kin.orientation();

			// Q2 and inelasticity cuts
			if (kin.Q2() > 4 * GeV2) vetoEvent;
			if (!inRange(kin.y(), 0.2, 0.8)) vetoEvent;

			// Jet calculation
			const Jets jets = apply<FastJets>(event, "Jets") \
				.jets(Cuts::Et > 6 * GeV && Cuts::etaIn(-1.375 * orientation, 1.875 * orientation), cmpMomByEt);
			MSG_DEBUG("Jet multiplicity = " << jets.size());
			//Dijet event selection
			if (jets.size() < 2) vetoEvent;
			const Jet & j1 = jets[0];
		    const Jet & j2 = jets[1];
	
			//Jet eta, average eta and eta difference calculation
			const double eta1 = orientation * j1.eta(), eta2 = orientation * j2.eta();
			const double etabar = (eta1 + eta2) / 2;
			const double etadiff = eta1 - eta2;

			//Cut in pseudorapidity difference
			if (abs(etadiff) > 0.5) vetoEvent;

			// Calculation of x_gamma^obs
			/// note Assuming Ee is the lab frame positron momentum, not in proton rest frame cf.
		    const double xyobs = (j1.Et() * exp(-eta1) + j2.Et() * exp(-eta2)) / (2 * kin.y() * kin.beamLepton().E());
			const size_t i_xyobs = (xyobs < 0.75) ? 0 : 1;

			//Classify events according to minimum jet energies
			size_t iE_min = 0;
			if (j1.Et() > 8 * GeV && j2.Et() > 8 * GeV)
				iE_min = 1;
			if (j1.Et() > 11 * GeV && j2.Et() > 11 * GeV)
				iE_min = 2;
			if (j1.Et() > 15 * GeV && j2.Et() > 15 * GeV)
				iE_min = 3;

			//Fill histograms
			for (size_t isel = 0; isel <= iE_min; ++isel) {

				//T1
				_h_etabar_all[isel]->fill(etabar);
				//T6, T7
				if (xyobs < 0.3) vetoEvent;
				_h_etabar[i_xyobs][isel]->fill(etabar);
			}
			

		}
		
		// Finalize
		void finalize() {
			const double sf = crossSection() / nanobarn / sumOfWeights();
			for (size_t ix = 0; ix < 2; ++ix) {
				for (auto& h : _h_etabar[ix]) scale(h, sf);

			}
			for (auto& h : _h_etabar_all) scale(h, sf);

		}


	private:
		//name Histograms
		Histo1DPtr _h_etabar_all[4], _h_etabar[2][4];

	};
	
	RIVET_DECLARE_PLUGIN(ZEUS_1997_I450085);

}

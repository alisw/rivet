// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

/// @brief Add a short analysis description here
class LHCF_2015_I1351909 : public Analysis {
public:

	/// Constructor
	DEFAULT_RIVET_ANALYSIS_CTOR(LHCF_2015_I1351909);

	static constexpr bool lhcf_like = true;
	static constexpr int ndecay = 1;
	static constexpr int nbeam = 2;
	static constexpr double D1_begin = 82000.; //mm 60000.; //mm
	static constexpr double D1_end = 82000; //mm 90000.; //mm
	static constexpr double IPtoLHCf = 141050.; //mm

	/// @name Analysis methods

	bool isParticleFromCollision(Particle p, vector<Particle> parents) {
		bool beam[nbeam]={false};

		if(parents.size()==nbeam) {
			for ( int ipar=0; ipar < nbeam; ++ipar )
				beam[ipar] = parents[ipar].genParticle()->is_beam();
			if(beam[0] && beam[1])
				return true;
		}

		return false;
	}

	bool isParticleFromDecay(Particle p, vector<Particle> parents) {
		if(parents.size()==ndecay)
			return true;
		else
			return false;
	}

	bool isDeviated(Particle p, Particle parent) { //Select/Remove particles decayed between IP and LHCf
		GenVertex* pv = p.genParticle()->production_vertex();
		assert(pv != NULL);

		const double decay_vertex = pv->position().z()/mm;

		const double parent_charge = PID::charge(parent.pid());
		const double descendant_charge = PID::charge(p.pid());

		if(parent_charge == 0) { //Particles produced by neutral parent decay
			if(descendant_charge == 0) {
				return false;
			} else {
				if(decay_vertex >= D1_end)
					return false;
				else
					return true; //Remove charged descendants produced from decay before end of D1
			}
		} else { //Particles produced by charged parent decay
			if(decay_vertex <= D1_begin) {
				if(descendant_charge == 0)
					return false;
				else
					return true; //Remove charged descendants produced from decay before end of D1
			} else {
				return true; //Remove particles produced by charged parent decay after begin of D1
			}
		}

		return false;
	}

	bool isSameParticle(Particle p1, Particle p2) {
		if(p1.pid() == p2.pid() &&
				mom(p1).t() == mom(p2).t() &&
				mom(p1).x() == mom(p2).x() &&
				mom(p1).y() == mom(p2).y() &&
				mom(p1).z() == mom(p2).z())
			return true;
		else
			return false;
	}

	bool isAlreadyProcessed(Particle p, vector<Particle> list) {
		for(unsigned int ipar=0; ipar<list.size(); ++ipar)
			if(isSameParticle(p, list[ipar]))
				return true;

		return false;
	}

	/// This method return a fake pseudorapidity to check id decayed particle is in LHCf acceptance
	double RecomputeEta(Particle p) {
		GenVertex* pv = p.genParticle()->production_vertex();

		const double x0 = pv->position().x()/mm;
		const double y0 = pv->position().y()/mm;
		const double z0 = pv->position().z()/mm;

		const double px = p.px()/MeV;
		const double py = p.py()/MeV;
		const double pz = abs(p.pz()/MeV);

		const double dist_to_lhcf = IPtoLHCf - z0;
		const double x1 = x0 + (dist_to_lhcf * px/pz);
		const double y1 = y0 + (dist_to_lhcf * py/pz);

		const double r = sqrt(pow(x1, 2.)+pow(y1, 2.));
		const double theta = atan(abs(r / IPtoLHCf));
		const double pseudorapidity = - log (tan (theta/2.) );

		return pseudorapidity;
	}

	/// Book histograms and initialise projections before the run
	void init() {

		// Initialise and register projections
		//      declare(FinalState("FS");
		addProjection(FinalState(), "FS");

		// Book histograms
		_h_n_en_eta1 = bookHisto1D(1, 1, 1);
		_h_n_en_eta2 = bookHisto1D(1, 1, 2);
		_h_n_en_eta3 = bookHisto1D(1, 1, 3);

	}

	/// Perform the per-event analysis
	void analyze(const Event& event) {

		const double weight = event.weight();

		const FinalState &fs = applyProjection<FinalState> (event, "FS");
		Particles fs_particles = fs.particles();

		vector<Particle> processed_parents;
		processed_parents.clear();

		for (Particle& p: fs_particles ) {

			if(p.pz()/GeV<0.) continue;

			double eta = 0.;
			double en = 0.;

			if(lhcf_like) {
				//======================================================================
				//========== LHCf-like analysis ========================================
				//======================================================================

				vector<Particle> parents = p.parents();

				if(isParticleFromCollision(p, parents)) { //Particles directly produced in collisions
					if(!PID::isHadron(p.pid())) continue; //Remove non-hadron particles
					if(PID::charge(p.pid()) != 0) continue; //Remove charged particles

					eta = p.eta();
					en = p.E()/GeV;
				} else if(isParticleFromDecay(p, parents)) { //Particles produced from decay
					GenVertex* pv = p.genParticle()->production_vertex();
					assert(pv != NULL);

					const double decay_vertex = pv->position().z()/mm;
					Particle parent = parents[0];

					if(decay_vertex < IPtoLHCf) { //If decay happens before LHCf we consider descendants
						if(!PID::isHadron(p.pid())) continue; //Remove non-hadron descendants
						if(isDeviated(p, parent)) continue; //Remove descendants deviated by D1

						eta = RecomputeEta(p);
						en = p.E()/GeV;
					} else {//If decay happens after LHCf we consider parents
						vector<Particle> ancestors;
						ancestors.clear();

						int ngeneration=0;
						bool isValid=true;
						bool isEnded=false;
						while(!isEnded) //Loop over all generations in the decay
						{
							vector<Particle> temp_part;
							temp_part.clear();
							if(ngeneration==0) {
								parent = parents[0];
								temp_part = parent.parents();
							}
							else {
								parent = ancestors[0];
								temp_part = parent.parents();
							}
							ancestors.clear();
							ancestors = temp_part;

							Particle ancestor = ancestors[0];

							if(isParticleFromCollision(parent, ancestors)) { //if we found first particles produced in collisions we consider them
								isEnded=true;

								if(!PID::isHadron(parent.pid())) isValid=false; //Remove non-hadron ancestors/parents
								if(PID::charge(parent.pid()) != 0) isValid=false; //Remove charged ancestors/parents
								if(isAlreadyProcessed(parent, processed_parents))
									isValid=false; //Remove already processed ancestors/parents when looping other descendants
								else
									processed_parents.push_back(parent); //Fill ancestors/parents in the list

								eta = parent.eta();
								en = parent.E()/GeV;
							} else if (isParticleFromDecay(parent, ancestors)) { //if we found first particles produced entering LHCf we consider them
								GenVertex* pv_prev = parent.genParticle()->production_vertex();
								assert(pv_prev != NULL);

								const double previous_decay_vertex = pv_prev->position().z()/mm;

								if(previous_decay_vertex < IPtoLHCf) {
									isEnded=true;

									if(!PID::isHadron(parent.pid())) isValid=false; //Remove non-hadron ancestors/parents
									if(isDeviated(parent, ancestor)) isValid=false; //Remove ancestors/parents deviated by D1
									if(isAlreadyProcessed(parent, processed_parents))
										isValid=false; //Remove already processed ancestors/parents when looping other descendants
									else
										processed_parents.push_back(parent); //Fill ancestors/parents in the list

									eta = RecomputeEta(parent);
									en = parent.E()/GeV;
								}
							} else { //This condition should never happen
								cout << "Looping over particles generation ended without match : Exit..." << endl;
								exit(EXIT_FAILURE);
							}

							++ngeneration;
						}

						if(!isValid) continue;
					}
				} else { //This condition should never happen
					cout << "Particle seems not to be produced in collision or decay : Exit..." << endl;
					exit(EXIT_FAILURE);
				}

			} else {
				//======================================================================
				//========== Only neutrons at IP =======================================
				//======================================================================

				vector<Particle> parents = p.parents();

				//if(isParticleFromCollision(p, parents)) { //Particles directly produced in collisions
					if(p.pid() != 2112 ) continue;

					eta = p.eta();
					en = p.E()/GeV;
				//}
			}

			// Fill histograms
			if( eta > 10.76 ){
				_h_n_en_eta1->fill( en , weight );

			}else if(eta > 8.99 && eta < 9.22){
				_h_n_en_eta2->fill( en , weight );

			}else if(eta > 8.81 && eta < 8.99){
				_h_n_en_eta3->fill( en , weight );

			}
		}
	}


	/// Normalise histograms etc., after the run
	void finalize() {

		scale(_h_n_en_eta1, crossSection()/millibarn/sumOfWeights()); // norm to cross section
		scale(_h_n_en_eta2, crossSection()/millibarn/sumOfWeights()); // norm to cross section
		scale(_h_n_en_eta3, crossSection()/millibarn/sumOfWeights()); // norm to cross section

	}

	//@}


private:


	/// @name Histograms
	//@{
	Histo1DPtr _h_n_en_eta1;
	Histo1DPtr _h_n_en_eta2;
	Histo1DPtr _h_n_en_eta3;
	//@}


};


// The hook for the plugin system
DECLARE_RIVET_PLUGIN(LHCF_2015_I1351909);


}

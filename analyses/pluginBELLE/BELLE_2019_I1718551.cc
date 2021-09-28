// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi, kaon and proton spectra at 10.58 GeV
  class BELLE_2019_I1718551 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2019_I1718551);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs, "FS");
      declare(Thrust(fs),"Thrust");
      declare(Beam(), "Beams");

      for(unsigned int ix=0;ix<6;++ix) {
	double xmin=0.05;
	for(unsigned int iy=0;iy<18;++iy) {
	  xmin+=0.05;
	  // pions
	  if(ix==0&&iy>15) continue;
	  Histo1DPtr temp1;
	  book(temp1,1,ix+1,iy+1);
	  _pion[ix].add(xmin,xmin+0.05,temp1);
	  // kaons
	  Histo1DPtr temp2;
	  book(temp2,2,ix+1,iy+1);
	  _kaon[ix].add(xmin,xmin+0.05,temp2);
	  if(iy>15 || (iy>14&&ix<3) || (iy>11&&ix==0)) continue;
	  // protons
	  Histo1DPtr temp3;
	  book(temp3,3,ix+1,iy+1);
	  _proton[ix].add(xmin,xmin+0.05,temp3);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      // get the thrust 
      const double tbins[6]={0.7,0.8,0.85,0.9,0.95,1.0};
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // find the thrust bin
      unsigned int ithrust=0;
      for(;ithrust<6;++ithrust)
	if(thrust.thrust()<=tbins[ithrust]) break;
      Vector3 a1 = thrust.thrustMajorAxis();
      Vector3 a2 = thrust.thrustMinorAxis();
      // loop over the charged hadrons
      const FinalState & fs = apply<FinalState>(event, "FS");
      for(const Particle & charged : fs.particles(Cuts::abspid==211 or
						  Cuts::abspid==321 or
						  Cuts::abspid==2212)) {
	double xE = charged.momentum().t()/meanBeamMom;
	double pT = sqrt(sqr(a1.dot(charged.momentum().p3()))+
			 sqr(a2.dot(charged.momentum().p3())));
	if(charged.abspid()==211)
	  _pion  [ithrust].fill(xE,pT);
	else if(charged.abspid()==321)
	  _kaon  [ithrust].fill(xE,pT);
	else if(charged.abspid()==2212)
	  _proton[ithrust].fill(xE,pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // scale factor including bin width in x_E
      double fact = crossSection()/femtobarn/sumOfWeights()/0.05;
      for(unsigned int ix=0;ix<6;++ix) {
	for(Histo1DPtr histo : _pion  [ix].histos()) scale(histo,fact);
	for(Histo1DPtr histo : _kaon  [ix].histos()) scale(histo,fact);
	for(Histo1DPtr histo : _proton[ix].histos()) scale(histo,fact);
      }
    }

    //@}


    /// @name Histograms
    //@{
    BinnedHistogram _pion[6],_kaon[6],_proton[6];
    //@}


  };


  DECLARE_RIVET_PLUGIN(BELLE_2019_I1718551);

}

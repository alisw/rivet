// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Upsilon(5S) -> Upsilon(nS) pi0 pi0
  class BELLE_2013_I1247463 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1247463);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==9000553);
      declare(ufs, "UFS");
      DecayedParticles UPS(ufs);
      UPS.addStable(PID::PI0);
      UPS.addStable(   553);
      UPS.addStable(100553);
      UPS.addStable(200553);
      declare(UPS, "UPS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the Upsilon(5S) decays
      static const map<PdgId,unsigned int> mode[3] = {{ {    553,1}, { 111,2}},
						      { { 100553,1}, { 111,2}},
						      { { 200553,1}, { 111,2}}};
      DecayedParticles UPS = apply<DecayedParticles>(event, "UPS");
      for(unsigned int ix=0;ix<UPS.decaying().size();++ix) {
	unsigned int imode;
	for(imode=0;imode<3;++imode) {
	  if(UPS.modeMatches(ix,3,mode[imode])) break;
	}
	if(imode>2) continue;
	const Particles & pi0 = UPS.decayProducts()[ix].at( 111);
	const Particle  & ups = UPS.decayProducts()[ix].at(553+imode*100000)[0];
	double mUpsPi[2];
	for(unsigned int iy=0;iy<2;++iy)
	  mUpsPi[iy] = (ups.momentum()+pi0[iy].momentum()).mass();
	if(mUpsPi[0]<mUpsPi[1]) swap(mUpsPi[0],mUpsPi[1]);
	_h[imode][0]->fill(mUpsPi[0]);
	_h[imode][1]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	_h[imode][2]->fill(mUpsPi[1]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	  normalize(_h[ix][iy]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1247463);

}

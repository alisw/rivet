// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> chi_c1 and chi_c2
  class BELLE_2016_I1408873 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2016_I1408873);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==300553 or Cuts::pid==9000553), "UPS");
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 ||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 20443);
      BB.addStable(   445);
      BB.addStable(   310);
      BB.addStable(   111);
      declare(BB, "BB");
      // book histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_spect[ix],1,1,1+ix);
	for(int iy=0;iy<5;++iy) {
	  book(_h_three[iy][ix],2,1+ix,1+iy);
	  if(iy<2)
	    book(_h_four[ix][iy],3,1+iy,1+ix);
	  else
	    book(_h_four[ix][iy],4,iy-1,1+ix);
	}
      }
      book(_wUps,"/TMP/Ups4");
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
        if (id == 20443 ||  id == 445 ) {
          unstable.push_back(p);
        }
        else if(!p.children().empty())
          findDecayProducts(p, unstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // upsilon spectrum
      for (const Particle& ups :  apply<UnstableParticles>(event, "UPS").particles()) {
	_wUps->fill();
	LorentzTransform cms_boost;
	if (ups.p3().mod() > 0.001)
	  cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	Particles unstable;
	// Find the decay products we want
	findDecayProducts(ups,unstable);
	for(const Particle & p : unstable) {
	  double modp = cms_boost.transform(p.momentum()).p3().mod();
	  if(p.pid()==20443) _h_spect[0]->fill(modp);
	  else               _h_spect[1]->fill(modp);
	}
      }
      // exclusive decay modes
      static const map<PdgId,unsigned int> & mode1   = { { 20443,1},{-211,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 20443,1},{ 211,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode2   = { {   445,1},{-211,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {   445,1},{ 211,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 20443,1},{ 211,1}, { 310,1}};
      static const map<PdgId,unsigned int> & mode3CC = { { 20443,1},{-211,1}, { 310,1}};
      static const map<PdgId,unsigned int> & mode4   = { {   445,1},{ 211,1}, { 310,1}};
      static const map<PdgId,unsigned int> & mode4CC = { {   445,1},{-211,1}, { 310,1}};
      static const map<PdgId,unsigned int> & mode5   = { { 20443,1},{ 111,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode5CC = { { 20443,1},{ 111,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode6   = { { 20443,1},{ 211,1}, {-211,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode6CC = { { 20443,1},{ 211,1}, {-211,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode7   = { {   445,1},{ 211,1}, {-211,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode7CC = { {   445,1},{ 211,1}, {-211,1}, {-321,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int sign = 1, imode =-1, ichi=0, iK=310,ipi=211;
	if(BB.decaying()[ix].abspid()==511) {
	  iK = 321;
	  ipi=-211;
	  if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
	    ichi=20443;
	    sign=1;
	    imode=0;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
	    ichi=20443;
	    sign=-1;
	    imode=0;
	  }
	  else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
	    ichi=445;
	    sign=1;
	    imode=1;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
	    ichi=445;
	    sign=-1;
	    imode=1;
	  }
	  else
	    continue;
	}
	else {
	  if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3)) {
	    ichi=20443;
	    sign=1;
	    imode=2;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3CC)) {
	    ichi=20443;
	    sign=-1;
	    imode=2;
	  }
	  else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode4)) {
	    ichi=445;
	    sign=1;
	    imode=3;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode4CC)) {
	    ichi=445;
	    sign=-1;
	    imode=3;
	  }
	  else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode5)) {
	    ichi=20443;
	    sign=1;
	    imode=4;
	    iK=321;
	    ipi=111;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode5CC)) {
	    ichi=20443;
	    sign=-1;
	    imode=4;
	    iK=321;
	    ipi=111;
	  }
	  else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,4,mode6)) {
	    ichi=20443;
	    sign=1;
	    imode=5;
	    iK=321;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,4,mode6CC)) {
	    ichi=20443;
	    sign=-1;
	    imode=5;
	    iK=321;
	  }
	  else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,4,mode7)) {
	    ichi=445;
	    sign=1;
	    imode=6;
	    iK=321;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,4,mode7CC)) {
	    ichi=445;
	    sign=-1;
	    imode=6;
	    iK=321;
	  }
	  else
	    continue;
	}
	if(ipi!=111) ipi*=sign;
	if(iK !=310) iK *=sign;
	const Particle & pip  = BB.decayProducts()[ix].at(ipi )[0];
	const Particle & K   = BB.decayProducts()[ix].at(iK  )[0];
	const Particle & chi = BB.decayProducts()[ix].at(ichi)[0];
	if(imode<5) {
	  _h_three[imode][0]->fill((pip.momentum()+K  .momentum()).mass());
	  _h_three[imode][1]->fill((chi.momentum()+pip.momentum()).mass());
	}
	else {
	  const Particle & pim  = BB.decayProducts()[ix].at(-ipi)[0];
	  FourMomentum ppi = pim.momentum()+pip.momentum();
	  imode -=5;
	  _h_four[imode][0]->fill((chi.momentum()+ppi).mass());
	  _h_four[imode][1]->fill((chi.momentum()+pim.momentum()).mass());
	  _h_four[imode][1]->fill((chi.momentum()+pip.momentum()).mass());
	  _h_four[imode][2]->fill((K.momentum()+ppi).mass());
	  _h_four[imode][3]->fill((K.momentum()+pim.momentum()).mass());
	  _h_four[imode][4]->fill(ppi.mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_spect[ix], 1e4/2./ *_wUps);
	for(unsigned int iy=0;iy<5;++iy) {
	  normalize(_h_three[iy][ix],1.,false);
	  normalize(_h_four [ix][iy],1.,false);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spect[2];
    Histo1DPtr _h_three[5][2];
    Histo1DPtr _h_four [2][5];
    CounterPtr _wUps;
    /// @}
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d03-x01-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d03-x02-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d03-x01-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d03-x02-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d04-x01-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d04-x02-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d04-x03-y01
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d04-x01-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d04-x02-y02
// BEGIN YODA_SCATTER2D_V2 /REF/BELLE_2016_I1408873/d04-x03-y02


  };


  RIVET_DECLARE_PLUGIN(BELLE_2016_I1408873);

}

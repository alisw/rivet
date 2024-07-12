// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Charm meson spectra in bottom decays
  class BABAR_2006_I719111 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I719111);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_Bm_D0      , 1, 1, 1);
      book(_h_Bm_Dbar0   , 1, 1, 2);
      book(_h_Bm_Dp      , 2, 1, 1);
      book(_h_Bm_Dm      , 2, 1, 2);
      book(_h_Bm_Dsp     , 3, 1, 1);
      book(_h_Bm_Dsm     , 3, 1, 2);
      book(_h_Bm_lam     , 4, 1, 1);
      book(_h_Bm_lbar    , 4, 1, 2);
      book(_h_Bbar0_D0   , 5, 1, 1);
      book(_h_Bbar0_Dbar0, 5, 1, 2);
      book(_h_Bbar0_Dp   , 6, 1, 1);
      book(_h_Bbar0_Dm   , 6, 1, 2);
      book(_h_Bbar0_Dsp  , 7, 1, 1);
      book(_h_Bbar0_Dsm  , 7, 1, 2);
      book(_h_Bbar0_lam  , 8, 1, 1);
      book(_h_Bbar0_lbar , 8, 1, 2);

      book(_c_Bm   ,"/TMP/Bm");
      book(_c_Bbar0,"/TMP/B0");
    }

    void findDecayProducts(Particle p,
			   Particles & D0, Particles & Dbar0,
			   Particles & Dp, Particles & Dm,
			   Particles & Dsp, Particles & Dsm,
			   Particles & lam, Particles & lbar) {
      for(const Particle & child : p.children()) {
	if(child.pid()==421)
	  D0.push_back(child);
	else if(child.pid()==-421)
	  Dbar0.push_back(child);
	else if(child.pid()==411)
	  Dp.push_back(child);
	else if(child.pid()==-411)
	  Dm.push_back(child);
	else if(child.pid()==431)
	  Dsp.push_back(child);
	else if(child.pid()==-431)
	  Dsm.push_back(child);
	else if(child.pid()==4122)
	  lam.push_back(child);
	else if(child.pid()==-4122)
	  lbar.push_back(child);
	else if(!child.children().empty())
	  findDecayProducts(child,D0,Dbar0,Dp,Dm, Dsp, Dsm, lam, lbar);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & p : ufs.particles(Cuts::abspid==511 || Cuts::abspid==521)) {
      	if(p.abspid()==511)
      	  _c_Bbar0->fill();
      	else if(p.abspid()==521)
      	  _c_Bm   ->fill();
      	Particles D0, Dbar0, Dp, Dm, Dsp, Dsm, lam, lbar;
      	findDecayProducts(p,D0,Dbar0,Dp,Dm, Dsp, Dsm, lam, lbar);
      	LorentzTransform boost =
      	  LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
      	if(p.pid()>0) {
      	  swap(D0, Dbar0);
      	  swap(Dp, Dm);
      	  swap(Dsp, Dsm);
      	  swap(lam, lbar);
      	}
      	for(const Particle & child : D0) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_D0   ->fill(pChild);
      	  else
      	    _h_Bbar0_D0->fill(pChild);
      	}
      	for(const Particle & child : Dbar0) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_Dbar0   ->fill(pChild);
      	  else
      	    _h_Bbar0_Dbar0->fill(pChild);
      	}
      	for(const Particle & child : Dp) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_Dp   ->fill(pChild);
      	  else
      	    _h_Bbar0_Dp->fill(pChild);
      	}
      	for(const Particle & child : Dm) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_Dm   ->fill(pChild);
      	  else
      	    _h_Bbar0_Dm->fill(pChild);
      	}
      	for(const Particle & child : Dsp) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_Dsp   ->fill(pChild);
      	  else
      	    _h_Bbar0_Dsp->fill(pChild);
      	}
      	for(const Particle & child : Dsm) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_Dsm   ->fill(pChild);
      	  else
      	    _h_Bbar0_Dsm->fill(pChild);
      	}
      	for(const Particle & child : lam) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_lam   ->fill(pChild);
      	  else
      	    _h_Bbar0_lam->fill(pChild);
      	}
      	for(const Particle & child : lbar) {
      	  double pChild = boost.transform(child.momentum()).p3().mod();
      	  if(p.abspid()==521)
      	    _h_Bm_lbar   ->fill(pChild);
      	  else
      	    _h_Bbar0_lbar->fill(pChild);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Bm_D0      ,100./ _c_Bm->val());
      scale(_h_Bm_Dbar0   ,100./ _c_Bm->val());
      scale(_h_Bm_Dp      ,100./ _c_Bm->val());
      scale(_h_Bm_Dm      ,100./ _c_Bm->val());
      scale(_h_Bm_Dsp     ,100./ _c_Bm->val());
      scale(_h_Bm_Dsm     ,100./ _c_Bm->val());
      scale(_h_Bm_lam     ,100./ _c_Bm->val());
      scale(_h_Bm_lbar    ,100./ _c_Bm->val());
      scale(_h_Bbar0_D0   ,100./ _c_Bbar0->val());
      scale(_h_Bbar0_Dbar0,100./ _c_Bbar0->val());
      scale(_h_Bbar0_Dp   ,100./ _c_Bbar0->val());
      scale(_h_Bbar0_Dm   ,100./ _c_Bbar0->val());
      scale(_h_Bbar0_Dsp  ,100./ _c_Bbar0->val());
      scale(_h_Bbar0_Dsm  ,100./ _c_Bbar0->val());
      scale(_h_Bbar0_lam  ,100./ _c_Bbar0->val());
      scale(_h_Bbar0_lbar ,100./ _c_Bbar0->val());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Bm_D0,_h_Bm_Dbar0,_h_Bm_Dp,_h_Bm_Dm,_h_Bm_Dsp,_h_Bm_Dsm,
      _h_Bm_lam,_h_Bm_lbar,_h_Bbar0_D0;
    Histo1DPtr _h_Bbar0_Dbar0,_h_Bbar0_Dp,_h_Bbar0_Dm,
      _h_Bbar0_Dsp,_h_Bbar0_Dsm,_h_Bbar0_lam,_h_Bbar0_lbar;

    CounterPtr _c_Bm ,_c_Bbar0;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2006_I719111);


}

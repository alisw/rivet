// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace {

unsigned int factorial(unsigned int n) {
  if(n<=1) return 0;
  unsigned int out=1;
  for (unsigned int i=n;i>0;--i) out *= i;
  return out;
}

double rapidityAxis(Rivet::Vector4 momentum, Rivet::Vector3 axis) {
  // assume axis is unit vector
  double plong = momentum.vector3().dot(axis);
  return 0.5*log((momentum.t()+plong)/(momentum.t()-plong));
}

double cosThetaStar(const Rivet::Particle & p1, const Rivet::Particle & p2, Rivet::Vector3 axis) {
  Rivet::FourMomentum psum = p1.momentum()+p2.momentum();
  Rivet::LorentzTransform trans(Rivet::LorentzTransform::mkObjTransformFromBeta(-psum.betaVec()));
  Rivet::FourMomentum m1 = trans.transform(p1.momentum());
  Rivet::FourVector a4;
  a4.setX(axis.x());
  a4.setY(axis.y());
  a4.setZ(axis.z());
  a4.setT(1.);
  a4 = trans.transform(a4);
  return fabs(m1.vector3().dot(a4.vector3())/m1.vector3().mod()/a4.vector3().mod());
}

}

namespace Rivet {


  /// @brief lambda anti-lambda correlations
  class OPAL_2000_I474010 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_2000_I474010);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs, "FS");
      FastJets durhamjets(fs, FastJets::DURHAM, 0.7);
      durhamjets.useInvisibles(JetAlg::Invisibles::ALL);
      declare(durhamjets, "DurhamJets");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_histLLBar    ,  1, 1, 1);
      book(_histLL       ,  1, 1, 2);
      book(_histLLBarCorr,  1, 1, 3);
      book(_histLLBar_2jet    ,  2, 1, 1);
      book(_histLL_2jet       ,  2, 1, 2);
      book(_histLLBarCorr_2jet,  2, 1, 3);
      book(_histLLBar_3jet    ,  3, 1, 1);
      book(_histLL_3jet       ,  3, 1, 2);
      book(_histLLBarCorr_3jet,  3, 1, 3);
      book(_histdcosTheta,  4, 1, 1);
      book(_histdy       ,  5, 1, 1);
      // weights
      book(_weight2Jet,"/TMP/W2Jet");
      book(_weight3Jet,"/TMP/W3Jet");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      const Thrust& thrust = applyProjection<Thrust>(event, "Thrust");
      Vector3 axis = thrust.thrustAxis();

      // Final state of unstable particles to get particle spectra
      const UnstableParticles & ufs = applyProjection<UnstableParticles>(event, "UFS");

      Particles lambda,lambdaBar;
      for (const Particle& p : ufs.particles()) {
        const int id = p.pid();
	if(id==3122)
	  lambda.push_back(p);
	else if(id==-3122)
	  lambdaBar.push_back(p);
      }

      // get the number of jets
      const FastJets& durjet = applyProjection<FastJets>(event, "DurhamJets");

      unsigned int njet = durjet.clusterSeq()->n_exclusive_jets_ycut(0.005);
      if(njet==2)
	_weight2Jet->fill();
      else if(njet==3)
	_weight3Jet->fill();

      // need a pair of baryons
      if(lambda.size()+lambdaBar.size()<2)
	vetoEvent;

      // multiplicities
      int nsame = (factorial(lambda.size())+factorial(lambdaBar.size()))/2;
      int ndiff = lambda.size()*lambdaBar.size();
      _histLLBar    ->fill(_histLLBar    ->bin(0).xMid(), ndiff       );
      _histLL       ->fill(_histLL       ->bin(0).xMid(),       nsame );
      _histLLBarCorr->fill(_histLLBarCorr->bin(0).xMid(),(ndiff-nsame));
      
      // multiplicities for different numbers of jets
      if(njet==2) {
	_histLLBar_2jet    ->fill(_histLLBar_2jet    ->bin(0).xMid(), ndiff       );
	_histLL_2jet       ->fill(_histLL_2jet       ->bin(0).xMid(),       nsame );
	_histLLBarCorr_2jet->fill(_histLLBarCorr_2jet->bin(0).xMid(),(ndiff-nsame));
      }
      else if(njet==3) {
	_histLLBar_3jet    ->fill(_histLLBar_3jet    ->bin(0).xMid(), ndiff       );
	_histLL_3jet       ->fill(_histLL_3jet       ->bin(0).xMid(),       nsame );
	_histLLBarCorr_3jet->fill(_histLLBarCorr_3jet->bin(0).xMid(),(ndiff-nsame));
      }
      // uncorrelated lambda
      for(unsigned int ix=0;ix<lambda.size();++ix) {
	double y1 = rapidityAxis(lambda[ix].momentum(),axis);
	for(unsigned int iy=ix+1;iy<lambda.size();++iy) {
	  double y2 = rapidityAxis(lambda[iy].momentum(),axis);
	  if(njet==2) _histdy->fill(fabs(y1-y2),-1.0);
	  double ctheta = cosThetaStar(lambda[ix],lambda[iy],axis);
	  _histdcosTheta->fill(ctheta,-1.0);
	}
      }
      // uncorrelated lambdabar
      for(unsigned int ix=0;ix<lambdaBar.size();++ix) {
	double y1 = rapidityAxis(lambdaBar[ix].momentum(),axis);
	for(unsigned int iy=ix+1;iy<lambdaBar.size();++iy) {
	  double y2 = rapidityAxis(lambdaBar[iy].momentum(),axis);
	  if(njet==2) _histdy->fill(fabs(y1-y2),-1.0);
	  double ctheta = cosThetaStar(lambdaBar[ix],lambdaBar[iy],axis);
	  _histdcosTheta->fill(ctheta,-1.0);
	}
      }
      // all lambda lambdabar
      for(unsigned int ix=0;ix<lambda.size();++ix) {
	double y1 = rapidityAxis(lambda[ix].momentum(),axis);
	for(unsigned int iy=0;iy<lambdaBar.size();++iy) {
	  double y2 = rapidityAxis(lambdaBar[iy].momentum(),axis);
	  if(njet==2) _histdy->fill(fabs(y1-y2));
	  double ctheta = cosThetaStar(lambda[ix],lambdaBar[iy],axis);
	  _histdcosTheta->fill(ctheta);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histLLBar,1./sumOfWeights());
      scale(_histLL,1./sumOfWeights());
      scale(_histLLBarCorr,1./sumOfWeights());
      scale(_histLLBar_2jet,1./ *_weight2Jet);
      scale(_histLL_2jet,1./ *_weight2Jet);
      scale(_histLLBarCorr_2jet,1./ *_weight2Jet);
      scale(_histLLBar_3jet,1./ *_weight3Jet);
      scale(_histLL_3jet,1./ *_weight3Jet);
      scale(_histLLBarCorr_3jet,1./ *_weight3Jet);
      normalize(_histdcosTheta, 1.);
      normalize(_histdy       , 1.);
    }

    //@}


    /// @name Histograms
    //@{

    CounterPtr _weight2Jet;
    CounterPtr _weight3Jet;

    Histo1DPtr _histLLBar;
    Histo1DPtr _histLL;
    Histo1DPtr _histLLBarCorr;

    Histo1DPtr _histLLBar_2jet;
    Histo1DPtr _histLL_2jet;
    Histo1DPtr _histLLBarCorr_2jet;

    Histo1DPtr _histLLBar_3jet;
    Histo1DPtr _histLL_3jet;
    Histo1DPtr _histLLBarCorr_3jet;

    Histo1DPtr _histdcosTheta;
    Histo1DPtr _histdy;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2000_I474010);


}

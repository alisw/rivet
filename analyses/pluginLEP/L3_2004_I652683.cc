// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"
#include "fastjet/EECambridgePlugin.hh"

namespace Rivet {


  /// Jet rates and event shapes at LEP I+II
  class L3_2004_I652683 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_2004_I652683);

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections to use
      const FinalState FS;
      declare(FS, "FS");
      declare(Beam(), "beams");
      const ChargedFinalState CFS;
      declare(CFS, "CFS");
      const Thrust thrust(FS);
      declare(thrust, "thrust");
      declare(ParisiTensor(FS), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");
      declare(InitialQuarks(), "initialquarks");
      FastJets jadeJets = FastJets(FS, FastJets::JADE, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      FastJets durhamJets = FastJets(FS, FastJets::DURHAM, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      FastJets cambridgeJets = FastJets(FS, FastJets::CAM, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jadeJets, "JadeJets");
      declare(durhamJets, "DurhamJets");

      // Book the histograms
      if (isCompatibleWithSqrtS(91.2)) {
        // z pole
        book(_h_Thrust_udsc             , 47, 1, 1);
        book(_h_Thrust_bottom           , 47, 1, 2);
        book(_h_heavyJetmass_udsc       , 48, 1, 1);
        book(_h_heavyJetmass_bottom     , 48, 1, 2);
        book(_h_totalJetbroad_udsc      , 49, 1, 1);
        book(_h_totalJetbroad_bottom    , 49, 1, 2);
        book(_h_wideJetbroad_udsc       , 50, 1, 1);
        book(_h_wideJetbroad_bottom     , 50, 1, 2);
        book(_h_Cparameter_udsc         , 51, 1, 1);
        book(_h_Cparameter_bottom       , 51, 1, 2);
        book(_h_Dparameter_udsc         , 52, 1, 1);
        book(_h_Dparameter_bottom       , 52, 1, 2);
        book(_h_Ncharged                , "/TMP/NCHARGED"     , 28, 1, 57);
        book(_h_Ncharged_udsc           , "/TMP/NCHARGED_UDSC", 28, 1, 57);
        book(_h_Ncharged_bottom         , "/TMP/NCHARGED_B"   , 27, 3, 57);
        book(_h_scaledMomentum          , 65, 1, 1);
        book(_h_scaledMomentum_udsc     , 65, 1, 2);
        book(_h_scaledMomentum_bottom   ,  65, 1, 3);
      }
      else if(sqrtS()/GeV < 90) {
        int i1(-1),i2(-1);
        if(isCompatibleWithSqrtS(41.4)) {
          i1=0; i2=1;
        }
        else if(isCompatibleWithSqrtS(55.3)) {
          i1=0; i2=2;
        }
        else if(isCompatibleWithSqrtS(65.4)) {
          i1=0; i2=3;
        }
        else if(isCompatibleWithSqrtS(75.7)) {
          i1=1; i2=1;
        }
        else if(isCompatibleWithSqrtS(82.3)) {
          i1=1; i2=2;
        }
        else if(isCompatibleWithSqrtS(85.1)) {
          i1=1; i2=3;
        }
        else
          MSG_ERROR("Beam energy not supported!");
        book(_h_thrust , 21+i1,1,i2);
        book(_h_rho    , 26+i1,1,i2);
        book(_h_B_T    , 31+i1,1,i2);
        book(_h_B_W    , 36+i1,1,i2);
      }
      else if(sqrtS()/GeV > 120) {
        int i1(-1),i2(-1);
        if(isCompatibleWithSqrtS(130.1)) {
          i1=0; i2=1;
        }
        else if(isCompatibleWithSqrtS(136.1)) {
          i1=0; i2=2;
        }
        else if(isCompatibleWithSqrtS(161.3)) {
          i1=0; i2=3;
        }
        else if(isCompatibleWithSqrtS(172.3)) {
          i1=1; i2=1;
        }
        else if(isCompatibleWithSqrtS(182.8)) {
          i1=1; i2=2;
        }
        else if(isCompatibleWithSqrtS(188.6)) {
          i1=1; i2=3;
        }
        else if(isCompatibleWithSqrtS(194.4)) {
          i1=2; i2=1;
        }
        else if(isCompatibleWithSqrtS(200.2)) {
          i1=2; i2=2;
        }
        else if(isCompatibleWithSqrtS(206.2)) {
          i1=2; i2=3;
        }
        else
          MSG_ERROR("Beam energy not supported!");
	book(_h_thrust , 23+i1,1,i2);
	book(_h_rho    , 28+i1,1,i2);
	book(_h_B_T    , 33+i1,1,i2);
	book(_h_B_W    , 38+i1,1,i2);
	book(_h_C      , 41+i1,1,i2);
	book(_h_D      , 44+i1,1,i2);
	book(_h_N      , "/TMP/NCHARGED", 22, 9, 53);
	book(_h_xi     , 66+i1,1,i2);
	// and the jets
	int i3 = 3*i1+i2;
	book(_h_y_2_JADE  ,   i3,1,1);
	book(_h_y_3_JADE  ,   i3,1,2);
	book(_h_y_4_JADE  ,   i3,1,3);
	book(_h_y_5_JADE  ,   i3,1,4);
	book(_h_y_2_Durham, 9+i3,1,1);
	book(_h_y_3_Durham, 9+i3,1,2);
	book(_h_y_4_Durham, 9+i3,1,3);
	book(_h_y_5_Durham, 9+i3,1,4);
	if(i3==8||i3==9) {
	  book(_h_y_2_Cambridge,11+i3,1,1);
	  book(_h_y_3_Cambridge,11+i3,1,2);
	  book(_h_y_4_Cambridge,11+i3,1,3);
	  book(_h_y_5_Cambridge,11+i3,1,4);
	}
      }

      book(_sumW_udsc, "_sumW_udsc");
      book(_sumW_b, "_sumW_b");
      book(_sumW_ch, "_sumW_ch");
      book(_sumW_ch_udsc, "_sumW_ch_udsc");
      book(_sumW_ch_b, "_sumW_ch_b");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beam average momentum
      const ParticlePair& beams = apply<Beam>(event, "beams").beams();
      const double beamMomentum = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;

      // InitialQuarks projection to have udsc events separated from b events
      /// @todo Yuck!!! Eliminate when possible...
      int iflav = 0;
      // only need the flavour at Z pole
      if(_h_Thrust_udsc) {
	int flavour = 0;
	const InitialQuarks& iqf = apply<InitialQuarks>(event, "initialquarks");
	Particles quarks;
	if ( iqf.particles().size() == 2 ) {
	  flavour = iqf.particles().front().abspid();
	  quarks  = iqf.particles();
	} else {
	  map<int, Particle> quarkmap;
	  for (const Particle& p : iqf.particles()) {
	    if (quarkmap.find(p.pid()) == quarkmap.end()) quarkmap[p.pid()] = p;
	    else if (quarkmap[p.pid()].E() < p.E()) quarkmap[p.pid()] = p;
	  }
	  double max_energy = 0.;
	  for (int i = 1; i <= 5; ++i) {
	    double energy = 0.;
	    if (quarkmap.find(i) != quarkmap.end())
	      energy += quarkmap[ i].E();
	    if (quarkmap.find(-i) != quarkmap.end())
	      energy += quarkmap[-i].E();
	    if (energy > max_energy)
	      flavour = i;
	  }
	  if (quarkmap.find(flavour) != quarkmap.end())
	    quarks.push_back(quarkmap[flavour]);
	  if (quarkmap.find(-flavour) != quarkmap.end())
	    quarks.push_back(quarkmap[-flavour]);
	}
	// Flavour label
	/// @todo Change to a bool?
	iflav = (flavour == PID::DQUARK || flavour == PID::UQUARK || flavour == PID::SQUARK || flavour == PID::CQUARK) ? 1 : (flavour == PID::BQUARK) ? 5 : 0;
      }
      // Update weight sums
      if (iflav == 1) {
        _sumW_udsc->fill();
      } else if (iflav == 5) {
        _sumW_b->fill();
      }
      _sumW_ch->fill();

      // Charged multiplicity
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      if(_h_Ncharged) _h_Ncharged->fill(cfs.size());
      if (iflav == 1) {
        _sumW_ch_udsc->fill();
        _h_Ncharged_udsc->fill(cfs.size());
      } else if (iflav == 5) {
        _sumW_ch_b->fill();
        _h_Ncharged_bottom->fill(cfs.size());
      }
      else if(_h_N) {
	_h_N->fill(cfs.size());
      }

      // Scaled momentum
      const Particles& chparticles = cfs.particlesByPt();
      for (const Particle& p : chparticles) {
        const Vector3 momentum3 = p.p3();
        const double mom = momentum3.mod();
        const double scaledMom = mom/beamMomentum;
        const double logScaledMom = std::log(scaledMom);
        if(_h_scaledMomentum) _h_scaledMomentum->fill(-logScaledMom);
        if (iflav == 1) {
          _h_scaledMomentum_udsc->fill(-logScaledMom);
        } else if (iflav == 5) {
          _h_scaledMomentum_bottom->fill(-logScaledMom);
        }
	else if(_h_xi) {
	  _h_xi->fill(-logScaledMom);
	}
      }

      // Thrust
      const Thrust& thrust = applyProjection<Thrust>(event, "thrust");
      if (iflav == 1) {
        _h_Thrust_udsc->fill(thrust.thrust());
      } else if (iflav == 5) {
        _h_Thrust_bottom->fill(thrust.thrust());
      }
      else if(_h_thrust) {
        _h_thrust->fill(1.-thrust.thrust());
      }

      // C and D Parisi parameters
      const ParisiTensor& parisi = applyProjection<ParisiTensor>(event, "Parisi");
      if (iflav == 1) {
        _h_Cparameter_udsc->fill(parisi.C());
        _h_Dparameter_udsc->fill(parisi.D());
      } else if (iflav == 5) {
        _h_Cparameter_bottom->fill(parisi.C());
        _h_Dparameter_bottom->fill(parisi.D());
      }
      else if(_h_C) {
	_h_C->fill(parisi.C());
	_h_D->fill(parisi.D());
      }

      // The hemisphere variables
      const Hemispheres& hemisphere = applyProjection<Hemispheres>(event, "Hemispheres");
      if (iflav == 1) {
        _h_heavyJetmass_udsc->fill(hemisphere.scaledM2high());
        _h_totalJetbroad_udsc->fill(hemisphere.Bsum());
        _h_wideJetbroad_udsc->fill(hemisphere.Bmax());
      } else if (iflav == 5) {
        _h_heavyJetmass_bottom->fill(hemisphere.scaledM2high());
        _h_totalJetbroad_bottom->fill(hemisphere.Bsum());
        _h_wideJetbroad_bottom->fill(hemisphere.Bmax());
      }
      else if (_h_rho) {
        _h_rho->fill(hemisphere.scaledM2high());
        _h_B_T->fill(hemisphere.Bsum());
        _h_B_W->fill(hemisphere.Bmax());
      }
      // jade jet rates
      if(_h_y_2_JADE) {
	const FastJets& jadejet = apply<FastJets>(event, "JadeJets");
	if (jadejet.clusterSeq()) {
	  const double y_23 = jadejet.clusterSeq()->exclusive_ymerge_max(2);
	  const double y_34 = jadejet.clusterSeq()->exclusive_ymerge_max(3);
	  const double y_45 = jadejet.clusterSeq()->exclusive_ymerge_max(4);
	  const double y_56 = jadejet.clusterSeq()->exclusive_ymerge_max(5);
	  for (size_t i = 0; i < _h_y_2_JADE->numBins(); ++i) {
	    double ycut = _h_y_2_JADE->bin(i).xMid();
	    double width = _h_y_2_JADE->bin(i).width();
	    if (y_23 < ycut) _h_y_2_JADE->fillBin(i,width);
	  }
	  for (size_t i = 0; i < _h_y_3_JADE->numBins(); ++i) {
	    double ycut = _h_y_3_JADE->bin(i).xMid();
	    double width = _h_y_3_JADE->bin(i).width();
	    if (y_34 < ycut && y_23 > ycut) {
	      _h_y_3_JADE->fillBin(i,width);
	    }
	  }
	  for (size_t i = 0; i < _h_y_4_JADE->numBins(); ++i) {
	    double ycut = _h_y_4_JADE->bin(i).xMid();
	    double width = _h_y_4_JADE->bin(i).width();
	    if (y_45 < ycut && y_34 > ycut) {
	      _h_y_4_JADE->fillBin(i,width);
	    }
	  }
	  for (size_t i = 0; i < _h_y_5_JADE->numBins(); ++i) {
	    double ycut = _h_y_5_JADE->bin(i).xMid();
	    double width = _h_y_5_JADE->bin(i).width();
	    if (y_56 < ycut && y_45 > ycut) {
	      _h_y_5_JADE->fillBin(i,width);
	    }
	  }
	}
      }
      // Durham jet rates
      if(_h_y_2_Durham) {
	const FastJets& durhamjet = apply<FastJets>(event, "DurhamJets");
	if (durhamjet.clusterSeq()) {
	  const double y_23 = durhamjet.clusterSeq()->exclusive_ymerge_max(2);
	  const double y_34 = durhamjet.clusterSeq()->exclusive_ymerge_max(3);
	  const double y_45 = durhamjet.clusterSeq()->exclusive_ymerge_max(4);
	  const double y_56 = durhamjet.clusterSeq()->exclusive_ymerge_max(5);
	  for (size_t i = 0; i < _h_y_2_Durham->numBins(); ++i) {
	    double ycut = _h_y_2_Durham->bin(i).xMid();
	    double width = _h_y_2_Durham->bin(i).width();
	    if (y_23 < ycut) _h_y_2_Durham->fillBin(i,width);
	  }
	  for (size_t i = 0; i < _h_y_3_Durham->numBins(); ++i) {
	    double ycut = _h_y_3_Durham->bin(i).xMid();
	    double width = _h_y_3_Durham->bin(i).width();
	    if (y_34 < ycut && y_23 > ycut) {
	      _h_y_3_Durham->fillBin(i,width);
	    }
	  }
	  for (size_t i = 0; i < _h_y_4_Durham->numBins(); ++i) {
	    double ycut = _h_y_4_Durham->bin(i).xMid();
	    double width = _h_y_4_Durham->bin(i).width();
	    if (y_45 < ycut && y_34 > ycut) {
	      _h_y_4_Durham->fillBin(i,width);
	    }
	  }
	  for (size_t i = 0; i < _h_y_5_Durham->numBins(); ++i) {
	    double ycut = _h_y_5_Durham->bin(i).xMid();
	    double width = _h_y_5_Durham->bin(i).width();
	    if (y_56 < ycut && y_45 > ycut) {
	      _h_y_5_Durham->fillBin(i,width);
	    }
	  }
	}
      }
      // Cambridge
      if(_h_y_2_Cambridge) {
	PseudoJets pjs;
	const FinalState& fs = applyProjection<FinalState>(event, "FS");
	for (size_t i = 0; i < fs.particles().size(); ++i) {
	  fastjet::PseudoJet pj = fs.particles()[i];
	  pjs.push_back(pj);
	}
	for (size_t i = 0; i < _h_y_2_Cambridge->numBins(); ++i) {
	  double ycut = _h_y_2_Cambridge->bin(i).xMid();
	  double width = _h_y_2_Cambridge->bin(i).width();
	  fastjet::EECambridgePlugin plugin(ycut);
	  fastjet::JetDefinition jdef(&plugin);
	  fastjet::ClusterSequence cseq(pjs, jdef);
	  unsigned int njet = cseq.inclusive_jets().size();
	  if(njet==2)
	    _h_y_2_Cambridge->fillBin(i,width);
	  else if(njet==3) {
	    if(i<_h_y_3_Cambridge->numBins()) _h_y_3_Cambridge->fillBin(i,width);
	  }
	  else if(njet==4) {
	    if(i<_h_y_4_Cambridge->numBins()) _h_y_4_Cambridge->fillBin(i,width);
	  }
	  else if(njet==5) {
	    if(i<_h_y_5_Cambridge->numBins()) _h_y_5_Cambridge->fillBin(i,width);
	  }
	}
      }
    }
    
    Scatter2DPtr convertHisto(unsigned int ix,unsigned int iy, unsigned int iz, Histo1DPtr histo) {
      Scatter2D temphisto(refData(ix, iy, iz));
      Scatter2DPtr mult;
      book(mult, ix, iy, iz);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	double y    = histo->bins()[b].area();
	double yerr = histo->bins()[b].areaErr();
	mult->addPoint(x, y, ex, make_pair(yerr,yerr));
      }
      return mult;
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Z pole plots
      if(_h_Thrust_udsc) {
        scale(_h_Thrust_udsc,  1/_sumW_udsc->sumW());
        scale(_h_heavyJetmass_udsc,  1/_sumW_udsc->sumW());
        scale(_h_totalJetbroad_udsc, 1/_sumW_udsc->sumW());
	scale(_h_wideJetbroad_udsc, 1/_sumW_udsc->sumW());
        scale(_h_Cparameter_udsc, 1/_sumW_udsc->sumW());
        scale(_h_Dparameter_udsc, 1/_sumW_udsc->sumW());
        scale(_h_Thrust_bottom, 1./_sumW_b->sumW());
        scale(_h_heavyJetmass_bottom, 1./_sumW_b->sumW());
        scale(_h_totalJetbroad_bottom, 1./_sumW_b->sumW());
	scale(_h_wideJetbroad_bottom, 1./_sumW_b->sumW());
        scale(_h_Cparameter_bottom, 1./_sumW_b->sumW());
        scale(_h_Dparameter_bottom, 1./_sumW_b->sumW());
        scale(_h_Ncharged, 1./_sumW_ch->sumW());
        scale(_h_Ncharged_udsc, 1./_sumW_ch_udsc->sumW());
        scale(_h_Ncharged_bottom, 1./_sumW_ch_b->sumW());
        convertHisto(59,1,1,_h_Ncharged       );
        convertHisto(59,1,2,_h_Ncharged_udsc  );
        convertHisto(59,1,3,_h_Ncharged_bottom);

        scale(_h_scaledMomentum, 1/_sumW_ch->sumW());
        scale(_h_scaledMomentum_udsc, 1/_sumW_ch_udsc->sumW());
        scale(_h_scaledMomentum_bottom, 1/_sumW_ch_b->sumW());
      }
      else {
	if(_h_thrust) normalize(_h_thrust);
	if(_h_rho) normalize(_h_rho);
	if(_h_B_T) normalize(_h_B_T);
	if(_h_B_W) normalize(_h_B_W);
	if(_h_C) normalize(_h_C);
	if(_h_D) normalize(_h_D);
	if(_h_N) normalize(_h_N);
	if(_h_xi) scale(_h_xi,1./sumOfWeights());
	
      
      Scatter2DPtr mult;
	if(_h_N) {
    if(isCompatibleWithSqrtS(130.1)) {
	    convertHisto(60, 1, 1, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(136.1)) {
	    convertHisto(60, 1, 2, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(161.3)) {
	    convertHisto(60, 1, 3, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(172.3)) {
	    convertHisto(61, 1, 1, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(182.8)) {
	    convertHisto(61, 1, 2, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(188.6)) {
	    convertHisto(61, 1, 3, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(194.4)) {
	    convertHisto(62, 1, 1, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(200.2)) {
	    convertHisto(62, 1, 2, _h_N);
	  }
	  else if(isCompatibleWithSqrtS(206.2)) {
	    convertHisto(62, 1, 3, _h_N);
	  }
	}
	// // the jets
	if(_h_y_2_JADE) {
	  scale(_h_y_2_JADE, 1./ sumOfWeights());
	  scale(_h_y_3_JADE, 1./ sumOfWeights());
	  scale(_h_y_4_JADE, 1./ sumOfWeights());
	  scale(_h_y_5_JADE, 1./ sumOfWeights());
	}
	if(_h_y_2_Durham) {
	  scale(_h_y_2_Durham, 1./ sumOfWeights());
	  scale(_h_y_3_Durham, 1./ sumOfWeights());
	  scale(_h_y_4_Durham, 1./ sumOfWeights());
	  scale(_h_y_5_Durham, 1./ sumOfWeights());
	}
	if(_h_y_2_Cambridge) {
	  scale(_h_y_2_Cambridge, 1./ sumOfWeights());
	  scale(_h_y_3_Cambridge, 1./ sumOfWeights());
	  scale(_h_y_4_Cambridge, 1./ sumOfWeights());
	  scale(_h_y_5_Cambridge, 1./ sumOfWeights());
	}
      }
    }


    /// Weight counters
    CounterPtr _sumW_udsc, _sumW_b, _sumW_ch, _sumW_ch_udsc, _sumW_ch_b;

    /// @name Histograms
    //@{
    // at the Z pole
    Histo1DPtr _h_Thrust_udsc, _h_Thrust_bottom;
    Histo1DPtr _h_heavyJetmass_udsc, _h_heavyJetmass_bottom;
    Histo1DPtr _h_totalJetbroad_udsc, _h_totalJetbroad_bottom;
    Histo1DPtr _h_wideJetbroad_udsc, _h_wideJetbroad_bottom;
    Histo1DPtr _h_Cparameter_udsc, _h_Cparameter_bottom;
    Histo1DPtr _h_Dparameter_udsc, _h_Dparameter_bottom;
    Histo1DPtr _h_Ncharged, _h_Ncharged_udsc, _h_Ncharged_bottom;
    Histo1DPtr _h_scaledMomentum, _h_scaledMomentum_udsc, _h_scaledMomentum_bottom;
    // at other enegies
    Histo1DPtr _h_thrust,_h_rho,_h_B_T,_h_B_W,_h_C,_h_D,_h_N,_h_xi;
    // and the jets
    Histo1DPtr _h_y_2_JADE,_h_y_3_JADE,_h_y_4_JADE,_h_y_5_JADE;
    Histo1DPtr _h_y_2_Durham,_h_y_3_Durham,_h_y_4_Durham,_h_y_5_Durham;
    Histo1DPtr _h_y_2_Cambridge,_h_y_3_Cambridge,_h_y_4_Cambridge,_h_y_5_Cambridge;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(L3_2004_I652683);

}

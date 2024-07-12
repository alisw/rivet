// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "fastjet/EECambridgePlugin.hh"

namespace Rivet {


  /// @brief event shapes vs thrust direction
  class DELPHI_2000_I522656 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_2000_I522656);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections.
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");
      declare(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
      declare(FastJets(fs, FastJets::JADE  , 0.7), "JadeJets"  );

      // book histograms
      vector<double> bins={0.00,0.12,0.24,0.36,0.48,0.60,0.72,0.84,0.96};
      // thrust angle binned
      unsigned int iy=1,ioff=0;
      for(unsigned int ix=0;ix<8;++ix) {
	{Histo1DPtr tmp; _h_EEC.add(bins[ix],bins[ix+1],book(tmp,21+ioff,1,iy));}
	{Histo1DPtr tmp;_h_AEEC.add(bins[ix],bins[ix+1],book(tmp,25+ioff,1,iy));}
	{Histo1DPtr tmp;_h_cone.add(bins[ix],bins[ix+1],book(tmp,29+ioff,1,iy));}
	if(ioff==0) {
	  Histo1DPtr tmp;_h_thrust.add(bins[iy],bins[iy+1],book(tmp,33,1,iy));
	}
	++iy;
	if(iy==3) {
	  ++ioff;
	  iy=1;
	}
      }
      // total values
      book(_h_EEC_all   , 3,1,1);
      book(_h_AEEC_all  , 4,1,1);
      book(_h_cone_all  , 5,1,1);
      book(_h_thrust_all, 6,1,1);
      book(_h_Oblateness, 7,1,1);
      book(_h_C         , 8,1,1);
      book(_h_heavy     , 9,1,1);
      book(_h_sum       ,10,1,1);
      book(_h_diff      ,11,1,1);
      book(_h_wide      ,12,1,1);
      book(_h_total     ,13,1,1);
      book(_h_jade      ,17,1,1);
      book(_h_dur       ,18,1,1);
      book(_h_cam       ,20,1,1);
      book(_h_bin,"/TMP/hbin",bins);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if ( fs.particles().size() < 2) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // angle bettwen thrust and beam
      double cosThrust = abs(beams.first.p3().unit().dot(thrust.thrustAxis()));
      _h_bin->fill(cosThrust);
      // thrust and related
      _h_thrust_all->fill(           1.-thrust.thrust());
      _h_thrust     .fill(cosThrust, 1.-thrust.thrust());
      _h_Oblateness->fill(thrust.oblateness() );
      
      // visible energy and make pseudojets
      double Evis = 0.0;
      PseudoJets pjs;
      for (const Particle& p : fs.particles()) {
        Evis += p.E();
	fastjet::PseudoJet pj = p;
	pjs.push_back(pj);
      }
      double Evis2 = sqr(Evis);
      // (A)EEC
      // Need iterators since second loop starts at current outer loop iterator, i.e. no "foreach" here!
      for (Particles::const_iterator p_i = fs.particles().begin(); p_i != fs.particles().end(); ++p_i) {
        for (Particles::const_iterator p_j = p_i; p_j != fs.particles().end(); ++p_j) {
	  if(p_i == p_j) continue; 
          const Vector3 mom3_i = p_i->momentum().p3();
          const Vector3 mom3_j = p_j->momentum().p3();
          const double energy_i = p_i->momentum().E();
          const double energy_j = p_j->momentum().E();
          const double thetaij = 180.*mom3_i.unit().angle(mom3_j.unit())/M_PI;
          double eec = (energy_i*energy_j) / Evis2;
	  eec *= 2.;
          _h_EEC_all->fill(           thetaij, eec);
          _h_EEC     .fill(cosThrust, thetaij, eec);
          if (thetaij <90.){
	    _h_AEEC_all->fill(           thetaij, -eec);
	    _h_AEEC     .fill(cosThrust, thetaij, -eec);
	  }
          else {
	    _h_AEEC_all->fill(          180.-thetaij, eec);
	    _h_AEEC     .fill(cosThrust,180.-thetaij, eec);
	  }
        }
      }
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      _h_heavy->fill(hemi.scaledM2high());
      _h_diff ->fill(hemi.scaledM2diff());
      _h_sum  ->fill(hemi.scaledM2low()+hemi.scaledM2high());
      _h_wide ->fill(hemi.Bmax() );
      _h_total->fill(hemi.Bsum() );
      // C-parameter
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _h_C->fill(parisi.C());
      // jets
      const FastJets&  durjet = apply<FastJets>(event, "DurhamJets");
      const FastJets& jadejet = apply<FastJets>(event, "JadeJets");
      if (durjet .clusterSeq()) _h_dur ->fill( durjet.clusterSeq()->exclusive_ymerge_max(2));
      if (jadejet.clusterSeq()) _h_jade->fill(jadejet.clusterSeq()->exclusive_ymerge_max(2));
      // Cambridge is more complicated, inclusive defn
      for (size_t i = 0; i < _h_cam->numBins(); ++i) {
	double ycut = _h_cam->bin(i).xMax();
	// 	  double width = _h_y_2_Cambridge->bin(i).width();
	fastjet::EECambridgePlugin plugin(ycut);
	fastjet::JetDefinition jdef(&plugin);
	fastjet::ClusterSequence cseq(pjs, jdef);
	unsigned int njet = cseq.inclusive_jets().size();
	if(njet==2) {
	  _h_cam->fill( _h_cam->bin(i).xMid());
	  break;
	}
      }
      // jet cone
      Vector3 jetAxis=thrust.thrustAxis();
      if(hemi.highMassDirection()) jetAxis *=-1.;
      for(const Particle & p : fs.particles()) {
	const double thetaij = 180.*jetAxis.angle(p.p3().unit())/M_PI;
	double jcef = p.E()/ Evis;
	_h_cone_all->fill(          thetaij,jcef);
	_h_cone     .fill(cosThrust,thetaij,jcef);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<8;++ix) {
	if(ix<2) scale(_h_thrust.histos()[ix],1./_h_bin->bins()[ix].area());
	scale(_h_EEC.histos()[ix],180./M_PI/_h_bin->bins()[ix].area());
	scale(_h_AEEC.histos()[ix],180./M_PI/_h_bin->bins()[ix].area());
	scale(_h_cone.histos()[ix],180./M_PI/_h_bin->bins()[ix].area());
      }

      
      // _h_thrust.scale(1./sumOfWeights(),this);
      // _h_EEC.scale(180./M_PI/sumOfWeights(),this);
      // _h_AEEC.scale(180./M_PI/sumOfWeights(),this);
      // _h_cone.scale(180./M_PI/sumOfWeights(),this);



      
      scale(_h_thrust_all, 1./sumOfWeights());
      scale(_h_EEC_all, 180./M_PI/sumOfWeights());
      scale(_h_AEEC_all, 180./M_PI/sumOfWeights());
      scale(_h_cone_all, 180./M_PI/sumOfWeights());
      scale(_h_Oblateness, 1./sumOfWeights());
      scale(_h_C         , 1./sumOfWeights());
      scale(_h_heavy     , 1./sumOfWeights());
      scale(_h_sum       , 1./sumOfWeights());
      scale(_h_diff      , 1./sumOfWeights());
      scale(_h_wide      , 1./sumOfWeights());
      scale(_h_total     , 1./sumOfWeights());
      scale(_h_dur       , 1./sumOfWeights());
      scale(_h_jade      , 1./sumOfWeights());
      scale(_h_cam       , 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_thrust_all,_h_EEC_all,_h_AEEC_all,_h_cone_all;
    Histo1DPtr _h_Oblateness,_h_C,_h_heavy,_h_sum,_h_diff,_h_wide,_h_total;
    Histo1DPtr _h_jade,_h_dur,_h_cam;
    BinnedHistogram _h_thrust,_h_EEC,_h_AEEC,_h_cone;
    Histo1DPtr _h_bin;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_2000_I522656);

}

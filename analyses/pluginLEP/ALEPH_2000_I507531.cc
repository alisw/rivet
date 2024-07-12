// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi, eta, eta', K0, lambda spectra
  class ALEPH_2000_I507531 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_2000_I507531);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      declare(FinalState()       ,  "FS");

      // Histograms
      // incl
      book(_h_pi0 ,  1,1,1);
      book(_h_eta ,  2,1,1);
      book(_h_etaP,  3,1,1);
      book(_h_K0  , 16,1,1);
      book(_h_lam , 17,1,1);
      // two jet
      book(_h_2_pi0 ,  4,1,1);
      book(_h_2_eta ,  5,1,1);
      book(_h_2_etaP,  6,1,1);
      book(_h_2_K0  , 18,1,1);
      book(_h_2_lam , 19,1,1);
      // three jet
      book(_h_3_pi0 [0],  7,1,1);
      book(_h_3_pi0 [1],  8,1,1);
      book(_h_3_pi0 [2],  9,1,1);
      book(_h_3_eta [0], 10,1,1);
      book(_h_3_eta [1], 11,1,1);
      book(_h_3_eta [2], 12,1,1);
      book(_h_3_etaP[0], 13,1,1);
      book(_h_3_etaP[1], 14,1,1);
      book(_h_3_etaP[2], 15,1,1);
      book(_h_3_K0  [0], 20,1,1);
      book(_h_3_K0  [1], 21,1,1);
      book(_h_3_K0  [2], 22,1,1);
      book(_h_3_lam [0], 23,1,1);
      book(_h_3_lam [1], 24,1,1);
      book(_h_3_lam [2], 25,1,1);
      book(_w2,"/TMP/W2");
      book(_w3,"/TMP/W3");
    }

    void findDecayProducts(const Particle & parent, Particles & decay) {
      for(const Particle & child : parent.children()) {
	if(child.children().empty()) {
	  decay.push_back(child);
	}
	else {
	  findDecayProducts(child,decay);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;

      Particles decay,fs;
      // unstable particles
      const UnstableParticles ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & part : ufs.particles(Cuts::pid==111 or Cuts::pid==221 or Cuts::pid==331 or
						Cuts::pid==310 or Cuts::abspid==3122)) {
	fs.push_back(part);
	findDecayProducts(part,decay);
      }
      // FS particles    
      for(const Particle & part : apply<FinalState>(event, "FS").particles()) {
	bool skip=false;
	for(const Particle & dec :decay) {
	  if(dec.genParticle()==part.genParticle()) {
	    skip=true;
	    break;
	  }
	}
	if(skip) continue;
	fs.push_back(part);
      }
      // Definition of the Durham algorithm
      fastjet::JetDefinition durham_def(fastjet::ee_kt_algorithm, fastjet::E_scheme, fastjet::Best);
      // pseudojets
      vector<fastjet::PseudoJet> input_particles;
      // Pseudo-jets from the non photons
      unsigned int ix=0;
      for (const Particle& p : fs ) {
        const FourMomentum p4 = p.momentum();
        input_particles.push_back(fastjet::PseudoJet(p4.px(), p4.py(), p4.pz(), p4.E()));
	input_particles.back().set_user_index(ix);
	++ix;
      }
      // cluster the jets
      fastjet::ClusterSequence clust_seq(input_particles, durham_def);
      PseudoJets jets = fastjet::sorted_by_E(clust_seq.exclusive_jets_ycut(0.01));
      if(jets.size()==2)      _w2->fill();
      else if(jets.size()==3) _w3->fill();
      ix=0;
      for(const Particle & part : fs) {
         double xE = part.momentum().E()/meanBeamMom;
         double xP = part.momentum().p3().mod()/meanBeamMom;
	 int ijet = jets.size()!=3 ? -1 : findJet(ix,jets);
	 if(part.pid()==111) {
	   _h_pi0->fill(xE);
	   if(jets.size()==2) {
	     _h_2_pi0->fill(xE);
	   }
	   else if(jets.size()==3) {
	     _h_3_pi0[ijet]->fill(xE);
	   }
	 }
	 else if(part.pid()==221) {
	   _h_eta->fill(xE);
	   if(jets.size()==2) {
	     _h_2_eta->fill(xE);
	   }
	   else if(jets.size()==3) {
	     _h_3_eta[ijet]->fill(xE);
	   }
	 }
	 else if(part.pid()==331) {
	   _h_etaP->fill(xE);
	   if(jets.size()==2) {
	     _h_2_etaP->fill(xE);
	   }
	   else if(jets.size()==3) {
	     _h_3_etaP[ijet]->fill(xE);
	   }
	 }
	 else if(part.pid()==310) {
	   double xi=-log(xP);
	   _h_K0->fill(xi);
	   if(jets.size()==2) {
	     _h_2_K0->fill(xi);
	   }
	   else if(jets.size()==3) {
	     _h_3_K0[ijet]->fill(xi);
	   }
	 }
	 else if(part.abspid()==3122) {
	   double xi=-log(xP);
	   _h_lam->fill(xi);
	   if(jets.size()==2) {
	     _h_2_lam->fill(xi);
	   }
	   else if(jets.size()==3) {
	     _h_3_lam[ijet]->fill(xi);
	   }
	 }
	 else {
	   break;
	 }
	 ix+=1;
      }
    }

    int findJet(int id, const PseudoJets & jets) {
      for(unsigned int ijet=0;ijet<jets.size();++ijet) {
	for(const PseudoJet & con : jets[ijet].constituents()) {
	  if(con.user_index()==id)
	    return ijet;
	}
      }
      return -1;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pi0  ,1./sumOfWeights());
      scale(_h_eta  ,1./sumOfWeights());
      scale(_h_etaP ,1./sumOfWeights());
      scale(_h_K0   ,1./sumOfWeights());
      scale(_h_lam  ,1./sumOfWeights());
      scale(_h_2_pi0  ,1./ *_w2);
      scale(_h_2_eta  ,1./ *_w2);
      scale(_h_2_etaP ,1./ *_w2);
      scale(_h_2_K0   ,1./ *_w2);
      scale(_h_2_lam  ,1./ *_w2);
      for(unsigned int ix=0;ix<3;++ix ) {
	scale(_h_3_pi0[ix] ,1./ *_w3);
	scale(_h_3_eta[ix] ,1./ *_w3);
	scale(_h_3_etaP[ix],1./ *_w3);
	scale(_h_3_K0[ix]  ,1./ *_w3);
	scale(_h_3_lam[ix] ,1./ *_w3);
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pi0     , _h_eta     , _h_etaP     , _h_K0     , _h_lam     ;
    Histo1DPtr _h_2_pi0   , _h_2_eta   , _h_2_etaP   , _h_2_K0   , _h_2_lam   ;
    Histo1DPtr _h_3_pi0[3], _h_3_eta[3], _h_3_etaP[3], _h_3_K0[3], _h_3_lam[3];
    CounterPtr _w2,_w3;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALEPH_2000_I507531);


}

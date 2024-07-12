// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D -> K pi pi dalitz decays
  class MARKIII_1987_I247266 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKIII_1987_I247266);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411 or
						Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      declare(DD, "DD");

      // histos
      book(_h_1_Kmpip,1,1,1);
      book(_h_1_pipi ,1,1,2);
      book(_h_1_Kmpi0,1,1,3);
      book(_dalitz1, "dalitz1",50,0.3,3.2,50,0.,2.);
      
      book(_h_2_K0pip,2,1,1);
      book(_h_2_pipi ,2,1,2);
      book(_h_2_K0pim,2,1,3);
      book(_dalitz2, "dalitz2",50,0.3,3.2,50,0.,2.);
      
      book(_h_3_K0pip,3,1,1);
      book(_h_3_pipi ,3,1,2);
      book(_h_3_K0pi0,3,1,3);
      book(_dalitz3, "dalitz3",50,0.3,3.2,50,0.,2.);
      
      book(_h_4_Kmpip[0],4,1,1);
      book(_h_4_pipi    ,4,1,2);
      book(_h_4_Kmpip[1],4,1,3);
      book(_dalitz4, "dalitz4",50,0.3,2.,50,0.,2.);
      
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 211,1}, {-321,1}, {111,1} };
      static const map<PdgId,unsigned int> & mode1CC = { {-211,1}, { 321,1}, {111,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 211,1}, {-211,1}, {310,1} };
      static const map<PdgId,unsigned int> & mode3   = { { 211,1}, { 111,1}, {310,1} };
      static const map<PdgId,unsigned int> & mode3CC = { {-211,1}, { 111,1}, {310,1} };
      static const map<PdgId,unsigned int> & mode4   = { { 211,2}, {-321,1} };
      static const map<PdgId,unsigned int> & mode4CC = { {-211,2}, { 321,1} };
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	int sign = DD.decaying()[ix].pid()/DD.decaying()[ix].abspid();
	if(DD.decaying()[ix].abspid()==421) {
	  if ( ( DD.decaying()[ix].pid()>0 && DD.modeMatches(ix,3,mode1  )) ||
	       ( DD.decaying()[ix].pid()<0 && DD.modeMatches(ix,3,mode1CC))) {
	    const Particle & pi0 = DD.decayProducts()[ix].at(      111)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    const Particle & Km  = DD.decayProducts()[ix].at(-sign*321)[0];
	    double mneut  = (Km.momentum()+pip.momentum()).mass2();
	    double mminus = (Km.momentum()+pi0.momentum()).mass2();
	    double mpipi  = (pip.momentum()+pi0.momentum()).mass2();
	    _h_1_Kmpip->fill(mneut );
	    _h_1_pipi ->fill(mpipi );
	    _h_1_Kmpi0->fill(mminus);
	    _dalitz1  ->fill(mminus,mpipi);
	  }
	  else if ( DD.modeMatches(ix,3,mode2  )) {
	    const Particle & K0  = DD.decayProducts()[ix].at(      310)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    const Particle & pim = DD.decayProducts()[ix].at(-sign*211)[0];
	    double mminus = (pim.momentum()+K0.momentum() ).mass2();
	    double mplus  = (pip.momentum()+K0.momentum() ).mass2();
	    double mpipi  = (pip.momentum()+pim.momentum()).mass2();
	    _h_2_K0pip->fill(mplus);
	    _h_2_K0pim->fill(mminus);
	    _h_2_pipi ->fill(mpipi);
	    _dalitz2  ->fill(mminus,mpipi); 
	  }
	}
	else if(DD.decaying()[ix].abspid()==411) {
	  if(DD.modeMatches(ix,3,mode3  ) || DD.modeMatches(ix,3,mode3CC)) {
	    const Particle & pi0 = DD.decayProducts()[ix].at(      111)[0];
	    const Particle & K0  = DD.decayProducts()[ix].at(      310)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    double mminus = (K0.momentum()+pip.momentum()).mass2();
	    double mplus  = (K0.momentum()+pi0.momentum()).mass2();
	    double mpipi  = (pip.momentum()+pi0.momentum()).mass2();
	    _h_3_K0pip->fill( mminus);
	    _h_3_pipi ->fill( mpipi );
	    _h_3_K0pi0->fill( mplus );
	    _dalitz3  ->fill(mplus,mpipi);
	  }
	  else if(DD.modeMatches(ix,3,mode4  ) || DD.modeMatches(ix,3,mode4CC)) {
	    const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	    const Particle  & Km  = DD.decayProducts()[ix].at(-sign*321)[0];
	    double mplus  = (Km.momentum() +pip[0].momentum()).mass2();
	    double mminus = (Km.momentum() +pip[1].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pip[1].momentum()).mass2();
	    if(mplus<mminus) swap(mplus,mminus);
	    _h_4_Kmpip[1]->fill(mminus);
	    _h_4_Kmpip[0]->fill(mplus );
	    _h_4_pipi    ->fill( mpipi);
	    _dalitz4     ->fill(mminus,mpipi);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_1_Kmpip);
      normalize(_h_1_pipi );
      normalize(_h_1_Kmpi0);
      normalize(_dalitz1  );
      
      normalize(_h_2_K0pip);
      normalize(_h_2_pipi );
      normalize(_h_2_K0pim);
      normalize(_dalitz2  );
      
      normalize(_h_3_K0pip);
      normalize(_h_3_pipi );
      normalize(_h_3_K0pi0);
      normalize(_dalitz3  );
      
      normalize(_h_4_Kmpip[0]);
      normalize(_h_4_pipi    );
      normalize(_h_4_Kmpip[1]);
      normalize(_dalitz4     );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_1_Kmpip, _h_1_pipi, _h_1_Kmpi0;
    Histo2DPtr _dalitz1;
    Histo1DPtr _h_2_K0pip, _h_2_pipi, _h_2_K0pim;
    Histo2DPtr _dalitz2;
    Histo1DPtr _h_3_K0pip, _h_3_pipi, _h_3_K0pi0;
    Histo2DPtr _dalitz3;
    Histo1DPtr _h_4_Kmpip[2], _h_4_pipi;
    Histo2DPtr _dalitz4;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKIII_1987_I247266);

}

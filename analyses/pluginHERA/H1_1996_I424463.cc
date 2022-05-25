// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Transverse momentum spectra of charged particles in DIS (H1 1996)
  class H1_1996_I424463 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(H1_1996_I424463);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      //declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");

      // Book histograms
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
      declare(ChargedFinalState(), "CFS");
      declare(FinalState(), "FS");

 
      // Book histograms
      // specify custom binning
      //book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      //book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      //book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      
      book(_NevAll, "TMP/Nev_all");
      _h_dndpt_high_eta_bin.resize(10);
      _h_dndpt_low_eta_bin.resize(10);
      _hdndeta_pt1_bin.resize(10);
      _hdndeta_bin.resize(10);
      _hdndptmax_low_eta_bin.resize(8);
      int ixx = 0 ;
      for (size_t ix = 0; ix < 10; ++ix) {
        book(_Nevt_after_cuts[ix], "TMP/Nevt_after_cuts" + to_string(ix));
        book(_h_dndpt_high_eta_bin[ix], ix+1, 1, 1);
        book(_h_dndpt_low_eta_bin[ix], ix+11, 1, 1);
        if (ix != 6 && ix != 9) {
          book(_hdndptmax_low_eta_bin[ixx], ixx+21, 1, 1);
          ixx=ixx+1;
        }
        book(_hdndeta_pt1_bin[ix], ix+29, 1, 1);
        book(_hdndeta_bin[ix],  ix+39, 1, 1);
      }


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");      
      const FinalState& fs = apply<FinalState>(event, "FS");
      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;
      double W2 = dk.W2()/GeV;

      // Momentum of the scattered lepton
      FourMomentum leptonMom = dl.out().momentum();
      double enel = leptonMom.E();
      double thel = 180.-leptonMom.angle(dl.in().momentum())/degree;

      _NevAll -> fill() ;
      
       // cout <<"enel/GeV = "<<enel/GeV<<", thel = "<<thel<<", y = "<<y<<", x = "<<x<<std::endl;
       bool cut = y > 0.05 && Q2 > 5. && Q2 < 100.&& enel > 12. && W2 > 4400. && thel > 157. && thel < 173.;
       if (!cut) vetoEvent;
      
       int ibin[10] ;
       for ( int i=0; i< 10; i++) {
       ibin[i] = 0; } 
      
       if(5.<Q2&&Q2<50.&&x>0.0001&&x<0.001)    ibin[0]=1; 
       if(5.<Q2&&Q2<10.&&x>0.0001&&x<0.0002)   ibin[1]=1;
       if(6.<Q2&&Q2<10.&&x>0.0002&&x<0.0005)   ibin[2]=1;
       if(10.<Q2&&Q2<20.&&x>0.0002&&x<0.0005)  ibin[3]=1;
       if(10.<Q2&&Q2<20.&&x>0.0005&&x<0.0008)  ibin[4]=1;
       if(10.<Q2&&Q2<20.&&x>0.0008&&x<0.0015)  ibin[5]=1;
       if(10.<Q2&&Q2<20.&&x>0.0015&&x<0.004)   ibin[6]=1;
       if(20.<Q2&&Q2<50.&&x>0.0005&&x<0.0014)  ibin[7]=1;
       if(20.<Q2&&Q2<50.&&x>0.0014&&x<0.003)   ibin[8]=1;
       if(20.<Q2&&Q2<50.&&x>0.003&&x<0.01)     ibin[9]=1;

       for ( int i=0; i< 10; i++) { if(ibin[i]==1) _Nevt_after_cuts[i] ->fill(); } 

      // Extract the particles other than the lepton
      Particles particles;
      particles.reserve(fs.particles().size());
      ConstGenParticlePtr dislepGP = dl.out().genParticle();
      for (const Particle& p : fs.particles()) {
        ConstGenParticlePtr loopGP = p.genParticle();
        if (loopGP == dislepGP) continue;
        particles.push_back(p);
      }

      // Boost to hadronic CM
      const LorentzTransform hcmboost = dk.boostHCM();
      
      int mult = 0 ;
      // Loop over the particles
      // long ncharged(0);
      double ptmax_high[10], ptmax_low[10] ;
      for ( int i=0; i< 10; i++) {
       ptmax_high[i] = 0.; ptmax_low[i] = 0.; } 
      double EtSum = 0;
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
         const Particle& p = particles[ip1];

         double eta = p.momentum().pseudorapidity();
         // Boost to hcm
         const FourMomentum hcmMom = hcmboost.transform(p.momentum());
   
         // apply safety cuts
         if (eta > -5  && eta < 10.){
            mult = mult + 1;              
	      double pThcm =hcmMom.pT() ; 
            double etahcm = hcmMom.pseudorapidity();
            // cout << " charge " << PID::charge(p.pid()) << endl;
            if (etahcm > 0. && etahcm < 2.0){ EtSum = EtSum + hcmMom.Et();}
            if (PID::charge(p.pid()) != 0) {
              if (etahcm > 0.5 && etahcm < 1.5){ 
               for ( int i=0; i< 10; i++) { 
                  if(ibin[i]==1) {
                    // cout << " fill histo low "<< i << " pt = " << pThcm << endl; 
                       _h_dndpt_low_eta_bin[i] ->fill(pThcm);
                       if(pThcm > ptmax_low[i] ) ptmax_low[i] = pThcm;
                    } 
               } 

              }   
              if (etahcm > 1.5 && etahcm < 2.5){
               for ( int i=0; i< 10; i++) { 
                  if(ibin[i]==1) 
                    _h_dndpt_high_eta_bin[i] ->fill(pThcm); 
                    if(pThcm > ptmax_high[i] ) ptmax_high[i] = pThcm;
               } 
              }                    
              for ( int i=0; i< 10; i++) { 
              if(ibin[i]==1) _hdndeta_bin[i] ->fill(etahcm); 
              if(ibin[i]==1 && pThcm > 1.) _hdndeta_pt1_bin[i] ->fill(etahcm); 
              }
           }  
         }  // end of loop over the particles
      }
      int ii=0;
      for ( int i=0; i< 10; i++) { 
         if (i != 6 && i != 9 ) { 
             if( ibin[i]==1 && EtSum > 6. ) {
              _hdndptmax_low_eta_bin[ii] ->fill(ptmax_low[i]); 
               // cout << " filling ptmax " << ii << "  " <<  i << endl; 
             }
            ii=ii+1;
         }
      }

    }
    /// Normalise histograms etc., after the run
    void finalize() {
    cout << " All events: " << _NevAll->val() << " after cuts: "<<  _Nevt_after_cuts[0]->val() << endl;
    cout << " cut1 events: " << _NevAll->val() << " after cuts: "<<  _Nevt_after_cuts[1]->val() << endl;
     

     int ii =0;
     for ( int i=0; i< 10; i++) { 
        if (_Nevt_after_cuts[i]->val()  != 0) {
            scale(_h_dndpt_high_eta_bin[i], 1./ *_Nevt_after_cuts[i]);  
            scale(_h_dndpt_low_eta_bin[i], 1./ *_Nevt_after_cuts[i]); 
            scale(_hdndeta_bin[i], 1./ *_Nevt_after_cuts[i]); 
            scale(_hdndeta_pt1_bin[i], 1./ *_Nevt_after_cuts[i]); }
        if ( i !=6 && i!=9) { 
           // if (_Nevt_after_cuts[i]->val()  != 0) scale(_hdndptmax_low_eta_bin[ii], 1./ *_Nevt_after_cuts[i]); ii=ii+1; }; 
            if (_Nevt_after_cuts[i]->val()  != 0) normalize(_hdndptmax_low_eta_bin[ii]); ii=ii+1; }; 
     } 
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    array<CounterPtr,10> _Nevt_after_cuts;
    vector<Histo1DPtr> _h_dndpt_low_eta_bin, _h_dndpt_high_eta_bin;
    vector<Histo1DPtr> _hdndeta_bin, _hdndeta_pt1_bin, _hdndptmax_low_eta_bin;
    CounterPtr _NevAll ;
   ///@}


  };


  DECLARE_RIVET_PLUGIN(H1_1996_I424463);

}

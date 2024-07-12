// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Charged particle production in deep-inelastic ep scattering at H1 
  class H1_2013_I1217865 : public Analysis {
  public:

    /// Constructor
    
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2013_I1217865);
    

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      //declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");

      // Book histograms
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
      declare(ChargedFinalState(), "CFS");
      declare(FinalState(), "FS");
      _h_dn_dpT_cen.resize(9);
      _h_dn_dpT_curr.resize(9);
      _h_dn_deta_soft.resize(9);
      _h_dn_deta_hard.resize(9);

      book(_h_dn_dpT_cen[0],19,1,1); 
      book(_h_dn_dpT_curr[0] ,20,1,1); 
      book(_h_dn_deta_soft[0],1,1,1);
      book(_h_dn_deta_hard[0],2,1,1);
      for (size_t ix = 0; ix < 9; ++ix) {
        book(_Nevt_after_cuts[ix], "TMP/Nevt_after_cuts" + to_string(ix));
        if (ix > 0 ) {
          book(_h_dn_dpT_cen[ix], ix+20, 1, 1);
	    book(_h_dn_dpT_curr[ix],ix+28,1,1);
	    book(_h_dn_deta_soft[ix],ix+2,1,1);
	    book(_h_dn_deta_hard[ix],ix+10,1,1);
        }
      }
   

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");      
      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;

      // Momentum of the scattered lepton
      FourMomentum leptonMom = dl.out().momentum();
      double enel = leptonMom.E();
      double thel = 180.-leptonMom.angle(dl.in().momentum())/degree;


      getLog()<<Log::DEBUG<<"enel/GeV = "<<enel/GeV<<", thel = "<<thel<<", y = "<<y<<", x = "<<x<<std::endl;
       bool cut = y > 0.05 && y < 0.6 && Q2 > 5. && Q2 < 100.;
       if (!cut) vetoEvent;
      

       int ibin[10] ;
       for ( int i=0; i< 9; i++) {
       ibin[i] = 0; } 

       ibin[0] = 1 ;
       if(5.<Q2&&Q2<10.&&x>0.0001&&x<0.00024)   ibin[1]=1;
       if(5.<Q2&&Q2<10.&&x>0.00024&&x<0.0005)   ibin[2]=1;
       if(5.<Q2&&Q2<10.&&x>0.0005&&x<0.002)     ibin[3]=1;
       if(10.<Q2&&Q2<20.&&x>0.0002&&x<0.00052)  ibin[4]=1;
       if(10.<Q2&&Q2<20.&&x>0.00052&&x<0.0011)  ibin[5]=1;
       if(10.<Q2&&Q2<20.&&x>0.0011&&x<0.0037)   ibin[6]=1;
       if(20.<Q2&&Q2<100.&&x>0.0004&&x<0.0017)  ibin[7]=1;
       if(20.<Q2&&Q2<100.&&x>0.0017&&x<0.01)    ibin[8]=1; 

       for ( int i=0; i< 9; i++) { if(ibin[i]==1) _Nevt_after_cuts[i] ->fill(); } 


      // Extract the particles other than the lepton
      Particles particles;
      particles.reserve(cfs.particles().size());
      ConstGenParticlePtr dislepGP = dl.out().genParticle();
      for (const Particle& p : cfs.particles()) {
        ConstGenParticlePtr loopGP = p.genParticle();
        if (loopGP == dislepGP) continue;
        particles.push_back(p);
      }

      // Boost to hadronic CM
      const LorentzTransform hcmboost = dk.boostHCM();
      
      int mult = 0 ;
      // Loop over the particles
      // long ncharged(0);
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
      const Particle& p = particles[ip1];

      double eta = p.momentum().pseudorapidity();
      double pT = p.momentum().pT()/GeV; 

      // Boost to hcm
      const FourMomentum hcmMom = hcmboost.transform(p.momentum());
   
      if (pT > 0.15 && eta > -2. && eta < 2.5){

        mult = mult + 1;
               
	  double pThcm =hcmMom.pT() ; 
        double etahcm = hcmMom.pseudorapidity();


        if (etahcm > 0. && etahcm < 1.5){ 

            _h_dn_dpT_cen[0]->fill(pThcm );
	      for ( int i=1; i< 9; i++) { 
                  if(ibin[i]==1) { _h_dn_dpT_cen[i]->fill(pThcm );}}
       }   
	   
       if (etahcm > 1.5 && etahcm < 5.){
           _h_dn_dpT_curr[0]->fill(pThcm ); 
	      for ( int i=1; i< 9; i++) { 
                  if(ibin[i]==1) { _h_dn_dpT_curr[i]->fill(pThcm );}}
       }
	   
       if(pThcm < 1.){
           _h_dn_deta_soft[0]->fill(etahcm );
	      for ( int i=1; i< 9; i++) { 
                  if(ibin[i]==1) {_h_dn_deta_soft[i]->fill(etahcm );}}
      }
 
      if(pThcm > 1. && pThcm < 10.){
            _h_dn_deta_hard[0]->fill(etahcm ); 
	      for ( int i=1; i< 9; i++) { 
                  if(ibin[i]==1) {_h_dn_deta_hard[i]->fill(etahcm );}}
 
       }       
    }  // if (etahcm > 0. && etahcm < 1.5){
  }  // end of loop over the particles



    }


    /// Normalise histograms etc., after the run
    void finalize() {
    //  normalize(_h_dn_dpT_cen);

     if (_Nevt_after_cuts[0]->val()  != 0)  scale(_h_dn_dpT_cen[0],  1./ *_Nevt_after_cuts[0]);
     if (_Nevt_after_cuts[0]->val()  != 0)  scale(_h_dn_dpT_curr[0],  1./ *_Nevt_after_cuts[0]);
     if (_Nevt_after_cuts[0]->val()  != 0)  scale(_h_dn_deta_soft[0],  1./ *_Nevt_after_cuts[0]);
     if (_Nevt_after_cuts[0]->val()  != 0)  scale(_h_dn_deta_hard[0],   1./ *_Nevt_after_cuts[0]); 
     
     
     for ( int i=1; i< 9; i++) { 
        if (_Nevt_after_cuts[i]->val()  != 0) {
           scale(_h_dn_dpT_cen[i],  1./ *_Nevt_after_cuts[i]);
           scale(_h_dn_dpT_curr[i], 1./ *_Nevt_after_cuts[i]);
           scale(_h_dn_deta_soft[i],1./ *_Nevt_after_cuts[i]);
           scale(_h_dn_deta_hard[i],1./ *_Nevt_after_cuts[i]);
        }
     }



    }
  private:

    /**
     *  Polar angle with right direction of the beam
     */
    inline double beamAngle(const FourVector& v, const bool & order) {
      double thel = v.polarAngle()/degree;
      if(thel<0.) thel+=180.;
      if(!order) thel = 180.-thel;
      return thel;
    }

    //@}


    /// @name Histograms
    //@{
      Histo1DPtr _h_dn_dpT_2r;
      Histo1DPtr _h_dn_dpT_2l;

      vector<Histo1DPtr> _h_dn_dpT_cen;
      vector<Histo1DPtr> _h_dn_dpT_curr;      
      vector<Histo1DPtr> _h_dn_deta_soft;      
      vector<Histo1DPtr> _h_dn_deta_hard;
      array<CounterPtr,9> _Nevt_after_cuts;
  

    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(H1_2013_I1217865);


}

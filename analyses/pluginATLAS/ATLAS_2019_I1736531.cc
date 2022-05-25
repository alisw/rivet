// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {       

  /// @brief Underlying event in Z events
  class ATLAS_2019_I1736531 : public Analysis {
  public:
  
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1736531); 

    void init() {
   
      // Get options from the new option system
      PdgId flav = (getOption("LMODE") == "EL")? PID::ELECTRON : PID::MUON;

      //Projections 
      FinalState fs;
      ZFinder zfinder(fs, Cuts::abseta<2.4 && Cuts::pT>25.0*GeV, flav, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder, "ZFinder");      
      ChargedFinalState cfs(zfinder.remainingFinalState() );
      declare(cfs, "cfs");
    
    //Histograms
    
       book(_p["pTsum_tow_zpt"] , 1, 1, 1);
       book(_p["pTsum_trv_zpt"] , 2, 1, 1);
       book(_p["pTsum_tmin_zpt"] , 3, 1, 1);
       book(_p["pTsum_tmax_zpt"] , 4, 1, 1);
       book(_p["pTsum_away_zpt"] , 5, 1, 1);
       book(_p["nch_tow_zpt"] , 6, 1, 1);
       book(_p["nch_trv_zpt"] , 7, 1, 1);
       book(_p["nch_tmin_zpt"] , 8, 1, 1);
       book(_p["nch_tmax_zpt"] , 9, 1, 1);
       book(_p["nch_away_zpt"] , 10, 1, 1);
       book(_p["pTmean_tow_zpt"] , 11, 1, 1);
       book(_p["pTmean_trv_zpt"] , 12, 1, 1);
       book(_p["pTmean_tmin_zpt"] , 13, 1, 1);
       book(_p["pTmean_tmax_zpt"] , 14, 1, 1);
       book(_p["pTmean_away_zpt"] , 15, 1, 1);
       book(_p["pTsum_tow_zpt_tlow"] , 16, 1, 1);
       book(_p["pTsum_trv_zpt_tlow"] , 17, 1, 1);
       book(_p["pTsum_tmin_zpt_tlow"] , 18, 1, 1);
       book(_p["pTsum_tmax_zpt_tlow"] , 19, 1, 1);
       book(_p["pTsum_away_zpt_tlow"] , 20, 1, 1);
       book(_p["nch_tow_zpt_tlow"] , 21, 1, 1);
       book(_p["nch_trv_zpt_tlow"] , 22, 1, 1);
       book(_p["nch_tmin_zpt_tlow"] , 23, 1, 1);
       book(_p["nch_tmax_zpt_tlow"] , 24, 1, 1);
       book(_p["nch_away_zpt_tlow"] , 25, 1, 1);
       book(_p["pTmean_tow_zpt_tlow"] , 26, 1, 1);
       book(_p["pTmean_trv_zpt_tlow"] , 27, 1, 1);
       book(_p["pTmean_tmin_zpt_tlow"] , 28, 1, 1);
       book(_p["pTmean_tmax_zpt_tlow"] , 29, 1, 1);
       book(_p["pTmean_away_zpt_tlow"] , 30, 1, 1);
       book(_p["pTsum_tow_zpt_thi"] , 31, 1, 1);
       book(_p["pTsum_trv_zpt_thi"] , 32, 1, 1);
       book(_p["pTsum_tmin_zpt_thi"] , 33, 1, 1);
       book(_p["pTsum_tmax_zpt_thi"] , 34, 1, 1);
       book(_p["pTsum_away_zpt_thi"] , 35, 1, 1);
       book(_p["nch_tow_zpt_thi"] , 36, 1, 1);
       book(_p["nch_trv_zpt_thi"] , 37, 1, 1);
       book(_p["nch_tmin_zpt_thi"] , 38, 1, 1);
       book(_p["nch_tmax_zpt_thi"] , 39, 1, 1);
       book(_p["nch_away_zpt_thi"] , 40, 1, 1);
       book(_p["pTmean_tow_zpt_thi"] , 41, 1, 1);
       book(_p["pTmean_trv_zpt_thi"] , 42, 1, 1);
       book(_p["pTmean_tmin_zpt_thi"] , 43, 1, 1);
       book(_p["pTmean_tmax_zpt_thi"] , 44, 1, 1);
       book(_p["pTmean_away_zpt_thi"] , 45, 1, 1);    
        
      for (size_t i_bin = 0; i_bin < 8; ++i_bin) {
            
        book(_h["pT[0]"+to_str(i_bin)] , 46 + i_bin, 1, 1); 
        book(_h["pT[1]"+to_str(i_bin)] , 54 + i_bin, 1, 1); 
        book(_h["pT[2]"+to_str(i_bin)] , 62 + i_bin, 1, 1); 
        book(_h["pT[3]"+to_str(i_bin)] , 70 + i_bin, 1, 1); 
        book(_h["pT[4]"+to_str(i_bin)] , 78 + i_bin, 1, 1);     
        book(_h["nch[0]"+to_str(i_bin)] , 86 + i_bin, 1, 1);        
        book(_h["nch[1]"+to_str(i_bin)] , 94 + i_bin, 1, 1); 
        book(_h["nch[2]"+to_str(i_bin)] , 102 + i_bin, 1, 1); 
        book(_h["nch[3]"+to_str(i_bin)] , 110 + i_bin, 1, 1); 
        book(_h["nch[4]"+to_str(i_bin)] , 118 + i_bin, 1, 1);    
        book(_h["pTsum[0]"+to_str(i_bin)] , 126 + i_bin, 1, 1);        
        book(_h["pTsum[1]"+to_str(i_bin)] , 134 + i_bin, 1, 1); 
        book(_h["pTsum[2]"+to_str(i_bin)] , 142 + i_bin, 1, 1); 
        book(_h["pTsum[3]"+to_str(i_bin)] , 150 + i_bin, 1, 1); 
        book(_h["pTsum[4]"+to_str(i_bin)] , 158 + i_bin, 1, 1); 
        book(_h["pTmean[0]"+to_str(i_bin)] , 166 + i_bin, 1, 1);        
        book(_h["pTmean[1]"+to_str(i_bin)] , 174 + i_bin, 1, 1); 
        book(_h["pTmean[2]"+to_str(i_bin)] , 182 + i_bin, 1, 1); 
        book(_h["pTmean[3]"+to_str(i_bin)] , 190 + i_bin, 1, 1); 
        book(_h["pTmean[4]"+to_str(i_bin)] , 198 + i_bin, 1, 1); 
        book(_h["pT_tlow[0]"+to_str(i_bin)] , 206 + i_bin, 1, 1);        
        book(_h["pT_tlow[1]"+to_str(i_bin)] , 214 + i_bin, 1, 1); 
        book(_h["pT_tlow[2]"+to_str(i_bin)] , 222 + i_bin, 1, 1); 
        book(_h["pT_tlow[3]"+to_str(i_bin)] , 230 + i_bin, 1, 1); 
        book(_h["pT_tlow[4]"+to_str(i_bin)] , 238 + i_bin, 1, 1);     
        book(_h["nch_tlow[0]"+to_str(i_bin)] , 246 + i_bin, 1, 1);        
        book(_h["nch_tlow[1]"+to_str(i_bin)] , 254 + i_bin, 1, 1); 
        book(_h["nch_tlow[2]"+to_str(i_bin)] , 262 + i_bin, 1, 1); 
        book(_h["nch_tlow[3]"+to_str(i_bin)] , 270 + i_bin, 1, 1); 
        book(_h["nch_tlow[4]"+to_str(i_bin)] , 278 + i_bin, 1, 1);    
        book(_h["pTsum_tlow[0]"+to_str(i_bin)] , 286 + i_bin, 1, 1);        
        book(_h["pTsum_tlow[1]"+to_str(i_bin)] , 294 + i_bin, 1, 1); 
        book(_h["pTsum_tlow[2]"+to_str(i_bin)] , 302 + i_bin, 1, 1); 
        book(_h["pTsum_tlow[3]"+to_str(i_bin)] , 310 + i_bin, 1, 1); 
        book(_h["pTsum_tlow[4]"+to_str(i_bin)] , 318 + i_bin, 1, 1); 
        book(_h["pTmean_tlow[0]"+to_str(i_bin)] , 326 + i_bin, 1, 1);        
        book(_h["pTmean_tlow[1]"+to_str(i_bin)] , 334 + i_bin, 1, 1); 
        book(_h["pTmean_tlow[2]"+to_str(i_bin)] , 342 + i_bin, 1, 1); 
        book(_h["pTmean_tlow[3]"+to_str(i_bin)] , 350 + i_bin, 1, 1); 
        book(_h["pTmean_tlow[4]"+to_str(i_bin)] , 358 + i_bin, 1, 1); 
        book(_h["pT_thi[0]"+to_str(i_bin)] , 366 + i_bin, 1, 1);        
        book(_h["pT_thi[1]"+to_str(i_bin)] , 374 + i_bin, 1, 1); 
        book(_h["pT_thi[2]"+to_str(i_bin)] , 382 + i_bin, 1, 1); 
        book(_h["pT_thi[3]"+to_str(i_bin)] , 390 + i_bin, 1, 1); 
        book(_h["pT_thi[4]"+to_str(i_bin)] , 398 + i_bin, 1, 1);     
        book(_h["nch_thi[0]"+to_str(i_bin)] , 406 + i_bin, 1, 1);        
        book(_h["nch_thi[1]"+to_str(i_bin)] , 414 + i_bin, 1, 1); 
        book(_h["nch_thi[2]"+to_str(i_bin)] , 422 + i_bin, 1, 1); 
        book(_h["nch_thi[3]"+to_str(i_bin)] , 430 + i_bin, 1, 1); 
        book(_h["nch_thi[4]"+to_str(i_bin)] , 438 + i_bin, 1, 1);    
        book(_h["pTsum_thi[0]"+to_str(i_bin)] , 446 + i_bin, 1, 1);        
        book(_h["pTsum_thi[1]"+to_str(i_bin)] , 454 + i_bin, 1, 1); 
        book(_h["pTsum_thi[2]"+to_str(i_bin)] , 462 + i_bin, 1, 1); 
        book(_h["pTsum_thi[3]"+to_str(i_bin)] , 470 + i_bin, 1, 1); 
        book(_h["pTsum_thi[4]"+to_str(i_bin)] , 478 + i_bin, 1, 1); 
        book(_h["pTmean_thi[0]"+to_str(i_bin)] , 486 + i_bin, 1, 1);        
        book(_h["pTmean_thi[1]"+to_str(i_bin)] , 494 + i_bin, 1, 1); 
        book(_h["pTmean_thi[2]"+to_str(i_bin)] , 502 + i_bin, 1, 1); 
        book(_h["pTmean_thi[3]"+to_str(i_bin)] , 510 + i_bin, 1, 1); 
        book(_h["pTmean_thi[4]"+to_str(i_bin)] , 518 + i_bin, 1, 1);                 
        }
 
    }

     
 // Perform the per-event analysis
    void analyze(const Event& event) {

      const double area = 5.*2./3.*M_PI;
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");

      if (zfinder.bosons().size() != 1) vetoEvent;
      double  Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
      double  Zphi  = zfinder.bosons()[0].momentum().phi();

     // Determine Zpt region histo to fill
     
      int i_bin(0);
      if (inRange(Zpt,0,10)) i_bin=0;
      if (inRange(Zpt,10,20)) i_bin=1;
      if (inRange(Zpt,20,40)) i_bin=2;
      if (inRange(Zpt,40,60)) i_bin=3;
      if (inRange(Zpt,60,80)) i_bin=4;
      if (inRange(Zpt,80,120)) i_bin=5;
      if (inRange(Zpt,120,200)) i_bin=6;
      if (Zpt>200) i_bin=7;


      // Initialization 
      int nTow(0), nTrans(0), nTransmin(0), nTransmax(0), nAway(0), nLeft(0), nRight(0);
      double pTsumTow(0.0), pTsumTrans(0.0), pTsumTransmin(0.0), pTsumTransmax(0.0), pTsumAway(0.0), pTsumLeft(0.0), pTsumRight(0.0);
    //  double pTmeanTow(0.0), pTmeanTrans(0.0), pTmeanTransmin(0.0), pTmeanTransmax(0.0), pTmeanAway(0.0);
      std::vector<double> leftpt;
      std::vector<double> rightpt;
  
      const Cut& pcut = ( (Cuts::abspid != PID::SIGMAMINUS) && (Cuts::abspid != PID::SIGMAPLUS) &&
                          (Cuts::abspid != PID::XIMINUS)    && (Cuts::abspid != PID::OMEGAMINUS) );

      Particles particles = apply<ChargedFinalState>(event, "cfs").particlesByPt(Cuts::pT > 0.5*GeV && Cuts::abseta <2.5 && pcut);


    //Calculate thrust
      vector<Vector3> momenta;
      for(const Particle& p : particles) {
          Vector3 mom = p.momentum().vector3();
           mom.setZ(0.0);
           momenta.push_back(mom);
           }

       if (momenta.size() == 2) {
          momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
        }
        
       Thrust thrustC;
       thrustC.calc(momenta);
       double thrust = thrustC.thrust();

      // Loop over charged particles 
      
        for(const Particle& p : particles) {
              double dphi = p.momentum().phi() - Zphi;
              double pT   = p.momentum().pT();
              for(; std::fabs(dphi) > M_PI; dphi += (dphi > 0. ? -2.*M_PI : 2.*M_PI) );

        // Towards region
        if( std::fabs(dphi) < M_PI/3. ) {
          nTow++;
          pTsumTow += pT;          
          _h["pT[0]"+to_str(i_bin)]->fill(pT);
          if(thrust < 0.75) _h["pT_tlow[0]"+to_str(i_bin)]->fill(pT);
          if(thrust > 0.75) _h["pT_thi[0]"+to_str(i_bin)]->fill(pT);    
        } 
        
        
        // Transverse region 
        else if( std::fabs(dphi) < 2.*M_PI/3. ) {        
          nTrans++;
          pTsumTrans += pT;          
          _h["pT[1]"+to_str(i_bin)]->fill(pT);
          if(thrust < 0.75) _h["pT_tlow[1]"+to_str(i_bin)]->fill(pT);
          if(thrust > 0.75) _h["pT_thi[1]"+to_str(i_bin)]->fill(pT);                    
           if(dphi > 0.) {
            nRight++;
            pTsumRight += pT;
            rightpt.push_back(pT);
          }
          else {
            nLeft++;
            pTsumLeft += pT;
            leftpt.push_back(pT);
          }
                                            
        }
        
                               
        // Away region
        else {
          nAway++;
          pTsumAway += pT;          
          _h["pT[4]"+to_str(i_bin)]->fill(pT);
          if(thrust < 0.75) _h["pT_tlow[4]"+to_str(i_bin)]->fill(pT);
          if(thrust > 0.75) _h["pT_thi[4]"+to_str(i_bin)]->fill(pT);          
        }        
     }
  
     // TransMAX, TransMIN regions
      if (pTsumLeft > pTsumRight) {
        pTsumTransmax = pTsumLeft;
        pTsumTransmin = pTsumRight;
        nTransmax     = nLeft;
        nTransmin     = nRight;
        
         for (auto x: rightpt) {       
             _h["pT[2]"+to_str(i_bin)]->fill(x);
            if(thrust < 0.75) _h["pT_tlow[2]"+to_str(i_bin)]->fill(x);
            if(thrust > 0.75) _h["pT_thi[2]"+to_str(i_bin)]->fill(x); 
             }        
         for (auto x: rightpt) {
            _h["pT[3]"+to_str(i_bin)]->fill(x);
            if(thrust < 0.75) _h["pT_tlow[3]"+to_str(i_bin)]->fill(x);
            if(thrust > 0.75) _h["pT_thi[3]"+to_str(i_bin)]->fill(x); 
             }       
       }
      
      else {
        pTsumTransmax = pTsumRight;
        pTsumTransmin = pTsumLeft;
        nTransmax     = nRight;
        nTransmin     = nLeft;
        
         for (auto x: leftpt) {
             _h["pT[2]"+to_str(i_bin)]->fill(x);
             if(thrust < 0.75) _h["pT_tlow[2]"+to_str(i_bin)]->fill(x);
             if(thrust > 0.75) _h["pT_thi[2]"+to_str(i_bin)]->fill(x); 
             }
        
         for (auto x: leftpt) {
             _h["pT[3]"+to_str(i_bin)]->fill(x);
            if(thrust < 0.75) _h["pT_tlow[3]"+to_str(i_bin)]->fill(x);
            if(thrust > 0.75) _h["pT_thi[3]"+to_str(i_bin)]->fill(x); 
             }       
      }

     // Fill rest of the histogtams

           
      _p["pTsum_tow_zpt"]->fill(Zpt, pTsumTow/area);
      _p["pTsum_trv_zpt"]->fill(Zpt, pTsumTrans/area);
      _p["pTsum_away_zpt"]->fill(Zpt, pTsumAway/area);
      _p["pTsum_tmin_zpt"]->fill(Zpt, pTsumTransmin/(0.5*area));
      _p["pTsum_tmax_zpt"]->fill(Zpt, pTsumTransmax/(0.5*area));     
      _p["nch_tow_zpt"]->fill(Zpt, nTow/area);
      _p["nch_trv_zpt"]->fill(Zpt, nTrans/area);
      _p["nch_away_zpt"]->fill(Zpt, nAway/area);
      _p["nch_tmin_zpt"]->fill(Zpt, nTransmin/(0.5*area));
      _p["nch_tmax_zpt"]->fill(Zpt, nTransmax/(0.5*area));
     
       if(nTow > 0)_p["pTmean_tow_zpt"]->fill(Zpt, pTsumTow/nTow);
       if(nTrans > 0)_p["pTmean_trv_zpt"]->fill(Zpt, pTsumTrans/nTrans);
       if(nAway > 0)_p["pTmean_away_zpt"]->fill(Zpt, pTsumAway/nAway);
       if(nTransmin > 0)_p["pTmean_tmin_zpt"]->fill(Zpt, pTsumTransmin/nTransmin);
       if(nTransmax > 0)_p["pTmean_tmax_zpt"]->fill(Zpt, pTsumTransmax/nTransmax);  
    
      if(thrust < 0.75){     
       _p["pTsum_tow_zpt_tlow"]->fill(Zpt, pTsumTow/area);
       _p["pTsum_trv_zpt_tlow"]->fill(Zpt, pTsumTrans/area);
       _p["pTsum_away_zpt_tlow"]->fill(Zpt, pTsumAway/area);
       _p["pTsum_tmin_zpt_tlow"]->fill(Zpt, pTsumTransmin/(0.5*area));
       _p["pTsum_tmax_zpt_tlow"]->fill(Zpt, pTsumTransmax/(0.5*area));     
       _p["nch_tow_zpt_tlow"]->fill(Zpt, nTow/area);
       _p["nch_trv_zpt_tlow"]->fill(Zpt, nTrans/area);
       _p["nch_away_zpt_tlow"]->fill(Zpt, nAway/area);
       _p["nch_tmin_zpt_tlow"]->fill(Zpt, nTransmin/(0.5*area));
       _p["nch_tmax_zpt_tlow"]->fill(Zpt, nTransmax/(0.5*area));

        if(nTow > 0)_p["pTmean_tow_zpt_tlow"]->fill(Zpt, pTsumTow/nTow);
        if(nTrans > 0)_p["pTmean_trv_zpt_tlow"]->fill(Zpt, pTsumTrans/nTrans);
        if(nAway > 0)_p["pTmean_away_zpt_tlow"]->fill(Zpt, pTsumAway/nAway);
        if(nTransmin > 0)_p["pTmean_tmin_zpt_tlow"]->fill(Zpt, pTsumTransmin/nTransmin);
        if(nTransmax > 0)_p["pTmean_tmax_zpt_tlow"]->fill(Zpt, pTsumTransmax/nTransmax);  
      }
     
     if(thrust > 0.75){          
        _p["pTsum_tow_zpt_thi"]->fill(Zpt, pTsumTow/area);
        _p["pTsum_trv_zpt_thi"]->fill(Zpt, pTsumTrans/area);
        _p["pTsum_away_zpt_thi"]->fill(Zpt, pTsumAway/area);
        _p["pTsum_tmin_zpt_thi"]->fill(Zpt, pTsumTransmin/(0.5*area));
        _p["pTsum_tmax_zpt_thi"]->fill(Zpt, pTsumTransmax/(0.5*area));
        _p["nch_tow_zpt_thi"]->fill(Zpt, nTow/area);
        _p["nch_trv_zpt_thi"]->fill(Zpt, nTrans/area);
        _p["nch_away_zpt_thi"]->fill(Zpt, nAway/area);
        _p["nch_tmin_zpt_thi"]->fill(Zpt, nTransmin/(0.5*area));
        _p["nch_tmax_zpt_thi"]->fill(Zpt, nTransmax/(0.5*area));

        if(nTow > 0)_p["pTmean_tow_zpt_thi"]->fill(Zpt, pTsumTow/nTow);
        if(nTrans > 0)_p["pTmean_trv_zpt_thi"]->fill(Zpt, pTsumTrans/nTrans);
        if(nAway > 0)_p["pTmean_away_zpt_thi"]->fill(Zpt, pTsumAway/nAway);
        if(nTransmin > 0)_p["pTmean_tmin_zpt_thi"]->fill(Zpt, pTsumTransmin/nTransmin);
        if(nTransmax > 0)_p["pTmean_tmax_zpt_thi"]->fill(Zpt, pTsumTransmax/nTransmax);     
           
      }
       
      _h["nch[0]"+to_str(i_bin)]->fill(nTow/area);
      _h["nch[1]"+to_str(i_bin)]->fill(nTrans/area);
      _h["nch[4]"+to_str(i_bin)]->fill(nAway/area);
      _h["nch[2]"+to_str(i_bin)]->fill(nTransmin/(0.5*area));
      _h["nch[3]"+to_str(i_bin)]->fill(nTransmax/(0.5*area));           
      _h["pTsum[0]"+to_str(i_bin)]->fill(pTsumTow/area);
      _h["pTsum[1]"+to_str(i_bin)]->fill(pTsumTrans/area);
      _h["pTsum[4]"+to_str(i_bin)]->fill(pTsumAway/area);
      _h["pTsum[2]"+to_str(i_bin)]->fill(pTsumTransmin/(0.5*area));
      _h["pTsum[3]"+to_str(i_bin)]->fill(pTsumTransmax/(0.5*area));

       if(nTow > 0) _h["pTmean[0]"+to_str(i_bin)]->fill(pTsumTow/nTow);
       if(nTrans > 0)_h["pTmean[1]"+to_str(i_bin)]->fill(pTsumTrans/nTrans);
       if(nAway > 0)_h["pTmean[4]"+to_str(i_bin)]->fill(pTsumAway/nAway);
       if(nTransmin > 0)_h["pTmean[2]"+to_str(i_bin)]->fill(pTsumTransmin/nTransmin);
       if(nTransmax > 0)_h["pTmean[3]"+to_str(i_bin)]->fill(pTsumTransmax/nTransmax);

       if(thrust < 0.75){      
       _h["nch_tlow[0]"+to_str(i_bin)]->fill(nTow/area);
       _h["nch_tlow[1]"+to_str(i_bin)]->fill(nTrans/area);
       _h["nch_tlow[4]"+to_str(i_bin)]->fill(nAway/area);
       _h["nch_tlow[2]"+to_str(i_bin)]->fill(nTransmin/(0.5*area));
       _h["nch_tlow[3]"+to_str(i_bin)]->fill(nTransmax/(0.5*area));    
       _h["pTsum_tlow[0]"+to_str(i_bin)]->fill(pTsumTow/area);
       _h["pTsum_tlow[1]"+to_str(i_bin)]->fill(pTsumTrans/area);
       _h["pTsum_tlow[4]"+to_str(i_bin)]->fill(pTsumAway/area);
       _h["pTsum_tlow[2]"+to_str(i_bin)]->fill(pTsumTransmin/(0.5*area));
       _h["pTsum_tlow[3]"+to_str(i_bin)]->fill(pTsumTransmax/(0.5*area));

       if(nTow > 0) _h["pTmean_tlow[0]"+to_str(i_bin)]->fill(pTsumTow/nTow);
       if(nTrans > 0)_h["pTmean_tlow[1]"+to_str(i_bin)]->fill(pTsumTrans/nTrans);
       if(nAway > 0)_h["pTmean_tlow[4]"+to_str(i_bin)]->fill(pTsumAway/nAway);
       if(nTransmin > 0)_h["pTmean_tlow[2]"+to_str(i_bin)]->fill(pTsumTransmin/nTransmin);
       if(nTransmax > 0)_h["pTmean_tlow[3]"+to_str(i_bin)]->fill(pTsumTransmax/nTransmax);
      }

      if(thrust > 0.75){
       _h["nch_thi[0]"+to_str(i_bin)]->fill(nTow/area);
       _h["nch_thi[1]"+to_str(i_bin)]->fill(nTrans/area);
       _h["nch_thi[4]"+to_str(i_bin)]->fill(nAway/area);
       _h["nch_thi[2]"+to_str(i_bin)]->fill(nTransmin/(0.5*area));
       _h["nch_thi[3]"+to_str(i_bin)]->fill(nTransmax/(0.5*area));            
       _h["pTsum_thi[0]"+to_str(i_bin)]->fill(pTsumTow/area);
       _h["pTsum_thi[1]"+to_str(i_bin)]->fill(pTsumTrans/area);
       _h["pTsum_thi[4]"+to_str(i_bin)]->fill(pTsumAway/area);
       _h["pTsum_thi[2]"+to_str(i_bin)]->fill(pTsumTransmin/(0.5*area));
       _h["pTsum_thi[3]"+to_str(i_bin)]->fill(pTsumTransmax/(0.5*area)); 

       if(nTow > 0) _h["pTmean_thi[0]"+to_str(i_bin)]->fill(pTsumTow/nTow);
       if(nTrans > 0)_h["pTmean_thi[1]"+to_str(i_bin)]->fill(pTsumTrans/nTrans);
       if(nAway > 0)_h["pTmean_thi[4]"+to_str(i_bin)]->fill(pTsumAway/nAway);
       if(nTransmin > 0)_h["pTmean_thi[2]"+to_str(i_bin)]->fill(pTsumTransmin/nTransmin);
       if(nTransmax > 0)_h["pTmean_thi[3]"+to_str(i_bin)]->fill(pTsumTransmax/nTransmax);
      }

  }
       
    
 void finalize() {
     
   normalize(_h);

 }

  private:
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
  };

  // This global object acts as a hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1736531);

}

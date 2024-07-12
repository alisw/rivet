// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/EventMixingFinalState.hh"
#include "Rivet/Projections/CentralityProjection.hh"
 #include "YODA/Utils/sortedvector.h"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class ALICE_2021_I1891391 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2021_I1891391);


    /// @name Analysis methods
    ///@{
    void fillbyparticles(BinnedHistogram &_histo, Particle tp, Particle ap) {
      double dphi = (ap.phi() - tp.phi());
      if (dphi < -0.5*M_PI) dphi = dphi + 2*M_PI;
      if (dphi > 1.5*M_PI) dphi = dphi - 2*M_PI;
      double deta = (ap.eta() - tp.eta());

      _histo.fill(deta, dphi, 1.0);
        
      return;
    }

    void S2DProjectionY(Scatter2DPtr projection, BinnedHistogram &hist2D) {

      double phiValues[N_phibins];
      double phiValueserr[N_phibins];

      vector<Point2D> points = projection->points();

      for (int i = 0; i < N_phibins; ++i)
      {
        phiValues[i]=0;
        phiValueserr[i]=0;
      }

      for(Histo1DPtr hist : hist2D.histos()) {
        int idx =0;
        for(auto bin:hist->bins()){
          phiValues[idx]+= bin.height();
          phiValueserr[idx]+=bin.heightErr()*bin.heightErr();
          idx+=1;
        }
      }
      projection->reset();

      for (int idx = 0; idx < N_phibins; ++idx) {
        phiValueserr[idx]=sqrt(phiValueserr[idx]);
        phiValues[idx]=phiValues[idx];
        projection->addPoint(points[idx].x(),phiValues[idx],points[idx].xErrAvg(),phiValueserr[idx]);
      }
      return;
    }

    pair<double, double> BackgEstimate (Scatter2DPtr hist) {

      pair<double, double> backg;
      vector<Point2D> points = hist->points();
      backg.first = (points[1].y()+points[2].y()+points[3].y()+points[34].y()+points[35].y()+points[36].y())/6;
      backg.second= pow(points[1].yErrAvg(),2)+pow(points[2].yErrAvg(),2)+pow(points[3].yErrAvg(),2);
      backg.second+= pow(points[34].yErrAvg(),2)+pow(points[35].yErrAvg(),2)+pow(points[36].yErrAvg(),2);

      backg.second=sqrt(backg.second)/6;

      return backg;
    }

    void ZYAM (Scatter2DPtr hist_final, Scatter2DPtr hist) {

      vector<Point2D> points = hist->points();
      
      pair<double, double> backg = BackgEstimate(hist);

      hist_final->reset();
      for(int idx = 0; idx < N_phibins; ++idx) {
        hist_final->addPoint(points[idx].x() ,points[idx].y()-backg.first,points[idx].xErrAvg(),sqrt(pow(points[idx].yErrAvg(),2)+pow(backg.second,2)));
      }
      return;
    }

    Point2D IntegratePeak (Scatter2DPtr s, pair<double, double> PeakInterval) {
      Point2D PeakYield;
      pair<double, double> Errs;
      PeakYield.setY(0);
      Errs.first = 0;
      Errs.second = 0;
      PeakYield.setYErrs(Errs);

      for (auto point : s->points()) {
        if (point.xMin() > PeakInterval.first && point.xMax() < PeakInterval.second) {
          PeakYield.setY(point.y() + PeakYield.y());
          Errs.first = sqrt(pow(point.yErrs().first,2) + pow(PeakYield.yErrs().first,2));
          Errs.second =  sqrt(pow(point.yErrs().second,2) + pow(PeakYield.yErrs().second,2));
          PeakYield.setYErrs(Errs);
        }
      }

      PeakYield.setY(PeakYield.y()*2*M_PI/N_phibins);
      Errs.first = PeakYield.yErrs().first*2*M_PI/N_phibins;
      Errs.second = PeakYield.yErrs().second*2*M_PI/N_phibins;
      PeakYield.setYErrs(Errs);

      return PeakYield;
    }

    void IntegratePeakByPT(Scatter2DPtr s, Scatter2DPtr vs[8], pair<double,double> PeakInterval, int pt_interval) {
      
      Point2D PeakYield;

      double xval_trig[PT_ASSOC_BINS] = { 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 8., 10., 13., 17.5};
      double xval_trigerr[PT_ASSOC_BINS] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1., 1., 2., 2.5};
      
      s->reset();
      int pt = 0;
      for (int pTsel = 0; pTsel < pt_interval+2; ++pTsel) {
        if(pt_interval==8){
            pt=pTsel+2;
            if(pTsel>7) continue; 
        }
        else pt=pTsel;
        PeakYield = IntegratePeak(vs[pTsel], PeakInterval);
        s->addPoint(xval_trig[pt],PeakYield.y(),xval_trigerr[pt],PeakYield.yErrAvg());
      }
      return;
    }

    void sdivide(Scatter2DPtr V0hist, Scatter2DPtr hhist, Scatter2DPtr ratio_hist,int pt_interval, bool V0mult) {
      vector<Point2D> V0h_points = V0hist->points();
      vector<Point2D> hh_points = hhist->points();

      double xval[PT_ASSOC_BINS] = { 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 8., 10., 13.,17.5};
      double xval_err[PT_ASSOC_BINS] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1., 1., 2., 2.5};
      double ratio,ratio_err;

      if(V0mult){
        xval[6]=9;
        xval[7]=15.5;
        xval_err[6]=2;
        xval_err[7]=4.5;
      }
      ratio_hist->reset();
      int pt = 0;
      for (int idx = 0; idx < pt_interval+2; ++idx) {
        if(pt_interval>=7||V0mult){
            pt = idx+2;
            if(idx>pt_interval-1) continue;
        }
        else pt = idx;
        if(V0mult&&idx>3){
          if(idx==4&&(hh_points[4].y()>0||hh_points[5].y())){
            ratio = (V0h_points[4].y()+V0h_points[5].y())/(hh_points[4].y()+hh_points[5].y());
            ratio_err=sqrt(pow(V0h_points[4].yErrAvg()/(hh_points[4].y()+hh_points[5].y()),2)+pow(V0h_points[5].yErrAvg()/(hh_points[4].y()+hh_points[5].y()),2)
                      +pow(((V0h_points[4].y()+V0h_points[5].y())*hh_points[4].yErrAvg())/(pow(hh_points[4].y()+hh_points[5].y(),2)),2)
                      +pow(((V0h_points[4].y()+V0h_points[5].y())*hh_points[5].yErrAvg())/(pow(hh_points[4].y()+hh_points[5].y(),2)),2));
          }
          else if(idx==5&&(hh_points[6].y()>0||hh_points[7].y())){
            ratio = (V0h_points[6].y()+V0h_points[7].y())/(hh_points[6].y()+hh_points[7].y());
            ratio_err=sqrt(pow(V0h_points[6].yErrAvg()/(hh_points[6].y()+hh_points[7].y()),2)+pow(V0h_points[7].yErrAvg()/(hh_points[6].y()+hh_points[7].y()),2)
                      +pow(((V0h_points[6].y()+V0h_points[7].y())*hh_points[6].yErrAvg())/(pow(hh_points[6].y()+hh_points[7].y(),2)),2)
                      +pow(((V0h_points[6].y()+V0h_points[7].y())*hh_points[7].yErrAvg())/(pow(hh_points[6].y()+hh_points[7].y(),2)),2));
          }else{
            ratio=1;
            ratio_err=0.02;
          }
        }
        else if(hh_points[idx].y()>0){
          ratio= V0h_points[idx].y()/hh_points[idx].y();
          ratio_err=sqrt(pow(V0h_points[idx].yErrAvg()/hh_points[idx].y(),2)+pow((V0h_points[idx].y()*hh_points[idx].yErrAvg())/(pow(hh_points[idx].y(),2)),2));
        }
        else {
          ratio=10;
          ratio_err=1;
        }
        ratio_hist->addPoint(xval[pt],ratio,xval_err[pt],ratio_err);
      }
    
      return;  
    }

    int profileIndex(vector<double> cBins, double c) {
      int index = 100;
      for (size_t i = 0; i < cBins.size() - 1; ++i) {
        if (c > cBins[i] && c <= cBins[i + 1]) {
          index = i;
          break;
        }
      }
      return index;
    }

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections

      //multiplicity
      declareCentrality(ALICE::V0MMultiplicity(),"ALICE_2015_PPCentrality","V0M","V0M");

      // Projections for trigger particles: charged, primary particles
      // with |eta| < 0.8 and different pT bins
      for (int ipt = 0; ipt < PT_TRIGG_BINS; ++ipt) {
        Cut cut = Cuts::abseta < 0.8 && Cuts::abscharge > 0 &&
          Cuts::ptIn(bins_pt_trigg[ipt]*GeV, bins_pt_trigg[ipt+1]*GeV);
        declare(ALICE::PrimaryParticles(cut), "APRIMTrigg" + toString(ipt));
      }

      // Projections for trigger particles: neutral, primary particles
      // with |y| < 0.5 and different pT bins
      for (int ipt = 0; ipt < PT_TRIGG_BINS; ++ipt) {
        Cut cut = Cuts::absrap < 0.5 && Cuts::abscharge == 0 &&
          Cuts::ptIn(bins_pt_trigg[ipt]*GeV, bins_pt_trigg[ipt+1]*GeV);
        declare(ALICE::PrimaryParticles(cut), "APRIMTrigg0" + toString(ipt));
      }

      // Projections for associated particles: charged, primary particles
      // with |eta| < 0.8 and different pT bins
      for (int ipt = 0; ipt < PT_ASSOC_BINS; ++ipt) {
        Cut cut = Cuts::abseta < 0.8 && Cuts::abscharge > 0 &&
          Cuts::ptIn(bins_pt_assoc[ipt]*GeV, bins_pt_assoc[ipt+1]*GeV);
        declare(ALICE::PrimaryParticles(cut), "APRIMAssoc" + toString(ipt));
      }
      const ChargedFinalState cfs(Cuts::pT > 1*GeV && Cuts::abseta < 0.8);
      declare(cfs, "CFS");
      const EventMixingFinalState evmc(cfs, cfs, 5, 0, 100, 10, 1.0);
      declare(evmc, "EVMc");

      multiplicityBins = {0.,1.,3.,7.,15.,50,100.};

      etabins[0]=-1.013333-2.666650e-02;
      for (int i = 1; i < Num_etabins; ++i){ etabins[i]=etabins[i-1]+(2*2.666650e-02);}
      
      // Histograms    

      for (int imult = 0; imult < MULT_BINS; ++imult) {
        for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg) {
          
          for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_hh_2D_mult[imult][ipt_trigg].add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_hh_2D_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg]+std::to_string(ieta), refData(2, 1, 1)));}
          for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_K0h_2D_mult[imult][ipt_trigg].add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_K0h_2D_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg]+std::to_string(ieta), refData(2, 1, 1)));}
          for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_Lamh_2D_mult[imult][ipt_trigg].add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_Lamh_2D_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg]+std::to_string(ieta), refData(2, 1, 1)));}
          
          book(_counterChargedTriggers_mult[imult][ipt_trigg], "TMP/counterChargedTriggers_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg]);
          book(_counterK0Triggers_mult[imult][ipt_trigg], "TMP/counterK0Triggers_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg]);
          book(_counterLamTriggers_mult[imult][ipt_trigg], "TMP/counterLamTriggers_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg]);
          
          book(_hist_dPhi_hh_mult[imult][ipt_trigg],"TMP/hist_dPhi_hh_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg], refData(2,1,1));
          book(_hist_dPhi_K0h_mult[imult][ipt_trigg],"TMP/hist_dPhi_K0h_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg], refData(3,1,1));
          book(_hist_dPhi_Lamh_mult[imult][ipt_trigg],"TMP/hist_dPhi_Lamh_mult_"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg], refData(4,1,1));

          if(imult<6||(imult==6&&!(ipt_trigg==0||ipt_trigg==5))){
            book(_hist_dPhi_hh_mult_fin[imult][ipt_trigg],"TMP/hist_dPhi_hh_mult_fin"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg], refData(2,1,1));
            book(_hist_dPhi_K0h_mult_fin[imult][ipt_trigg],"TMP/hist_dPhi_K0h_mult_fin"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg], refData(3,1,1));
            book(_hist_dPhi_Lamh_mult_fin[imult][ipt_trigg],"TMP/hist_dPhi_Lamh_mult_fin"+mulsel_name[imult]+"_"+bins_pt_trigg_name[ipt_trigg], refData(4,1,1));
          }
        } 
      }
      book(_hist_dPhi_hh_mult_fin[6][0],2,1,1);
      book(_hist_dPhi_K0h_mult_fin[6][0],3,1,1);
      book(_hist_dPhi_Lamh_mult_fin[6][0],4,1,1);

      book(_hist_dPhi_hh_mult_fin[6][5],5,1,1);
      book(_hist_dPhi_K0h_mult_fin[6][5],6,1,1);
      book(_hist_dPhi_Lamh_mult_fin[6][5],7,1,1);
         
      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg) {

        for (int ipt_assoc = 0; ipt_assoc < PT_ASSOC_BINS; ++ipt_assoc) {  
          for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_hh_2D_ptassoc[ipt_trigg][ipt_assoc].add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_hh_2D_ptassoc_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc]+std::to_string(ieta), refData(2, 1, 1)));}
          for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_K0h_2D_ptassoc[ipt_trigg][ipt_assoc].add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_K0h_2D_ptassoc_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc]+std::to_string(ieta), refData(2, 1, 1)));}
          for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_Lamh_2D_ptassoc[ipt_trigg][ipt_assoc].add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_Lamh_2D_ptassoc_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc]+std::to_string(ieta), refData(2, 1, 1)));}
          
          book(_hist_dPhi_hh_ptassoc[ipt_trigg][ipt_assoc],"TMP/hist_dPhi_hh_ptassoc_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc], refData(2,1,1));
          book(_hist_dPhi_K0h_ptassoc[ipt_trigg][ipt_assoc],"TMP/hist_dPhi_K0h_ptassoc_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc], refData(3,1,1));
          book(_hist_dPhi_Lamh_ptassoc[ipt_trigg][ipt_assoc],"TMP/hist_dPhi_Lamh_ptassoc_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc], refData(4,1,1));

          book(_hist_dPhi_hh_ptassoc_fin[ipt_trigg][ipt_assoc],"TMP/hist_dPhi_hh_ptassoc_fin_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc], refData(2,1,1));
          book(_hist_dPhi_K0h_ptassoc_fin[ipt_trigg][ipt_assoc],"TMP/hist_dPhi_K0h_ptassoc_fin_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc], refData(3,1,1));
          book(_hist_dPhi_Lamh_ptassoc_fin[ipt_trigg][ipt_assoc],"TMP/hist_dPhi_Lamh_ptassoc_fin_"+bins_pt_trigg_name[ipt_trigg]+bins_pt_assoc_name[ipt_assoc], refData(4,1,1));
          
        }
      }

      for (int imult = 0; imult < MULT_BINS; ++imult)
      {
        book(_hist_hh_NearSideYield_mult[imult],8,1,imult+1);
        book(_hist_K0h_NearSideYield_mult[imult],9,1,imult+1);
        book(_hist_Lamh_NearSideYield_mult[imult],10,1,imult+1);

        book(_hist_hh_AwaySideYield_mult[imult],11,1,imult+1);
        book(_hist_K0h_AwaySideYield_mult[imult],12,1,imult+1);
        book(_hist_Lamh_AwaySideYield_mult[imult],13,1,imult+1);

        

        if(imult<MULT_BINS-1){
          book(_hist_hh_NearSideYield_mult_ratio[imult],14,1,imult+1);
          book(_hist_K0h_NearSideYield_mult_ratio[imult],15,1,imult+1);
          book(_hist_Lamh_NearSideYield_mult_ratio[imult],16,1,imult+1);

          book(_hist_hh_AwaySideYield_mult_ratio[imult],17,1,imult+1);
          book(_hist_K0h_AwaySideYield_mult_ratio[imult],18,1,imult+1);
          book(_hist_Lamh_AwaySideYield_mult_ratio[imult],19,1,imult+1);

          book(_hist_K0h_NearSideYield_ratioTo_hh_mult[imult],imult+27,1,1);
          book(_hist_Lamh_NearSideYield_ratioTo_hh_mult[imult],imult+34,1,1);
          book(_hist_K0h_AwaySideYield_ratioTo_hh_mult[imult],imult+41,1,1);
          book(_hist_Lamh_AwaySideYield_ratioTo_hh_mult[imult],imult+48,1,1);
        }else{
          book(_hist_K0h_NearSideYield_ratioTo_hh_mult[imult],26,1,1);
          book(_hist_Lamh_NearSideYield_ratioTo_hh_mult[imult],33,1,1);
          book(_hist_K0h_AwaySideYield_ratioTo_hh_mult[imult],40,1,1);
          book(_hist_Lamh_AwaySideYield_ratioTo_hh_mult[imult],47,1,1);
        } 
      }
      
      for (int ipttrigg = 0; ipttrigg < PT_TRIGG_BINS; ++ipttrigg)
      {
        book(_hist_hh_NearSideYield_ptassoc[ipttrigg],20,1,ipttrigg+1);
        book(_hist_hh_AwaySideYield_ptassoc[ipttrigg],23,1,ipttrigg+1);

        book(_hist_K0h_NearSideYield_ptassoc[ipttrigg],21,1,ipttrigg+1);
        book(_hist_K0h_AwaySideYield_ptassoc[ipttrigg],24,1,ipttrigg+1);



        if(ipttrigg<PT_TRIGG_BINS-1){
          book(_hist_Lamh_NearSideYield_ptassoc[ipttrigg],22,1,ipttrigg+1);
          book(_hist_Lamh_AwaySideYield_ptassoc[ipttrigg],25,1,ipttrigg+1);

          book(_hist_K0h_NearSideYield_ratioTo_hh_ptassoc[ipttrigg],ipttrigg+54,1,1);
          book(_hist_Lamh_NearSideYield_ratioTo_hh_ptassoc[ipttrigg],ipttrigg+61,1,1);
          book(_hist_K0h_AwaySideYield_ratioTo_hh_ptassoc[ipttrigg],ipttrigg+68,1,1);
          book(_hist_Lamh_AwaySideYield_ratioTo_hh_ptassoc[ipttrigg],ipttrigg+75,1,1);
        }
      }

      for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_mix_hh.add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_mix_hh"+std::to_string(ieta), refData(2, 1, 1)));}
      for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_mix_K0h.add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_mix_K0h"+std::to_string(ieta), refData(2, 1, 1)));}
      for(int ieta=1; ieta<Num_etabins; ieta++){Histo1DPtr tmp; _hist_mix_Lamh.add(etabins[ieta-1], etabins[ieta], book(tmp, "TMP/hist_mix_Lamh"+std::to_string(ieta), refData(2, 1, 1)));}
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) { 

      //multiplicity block
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0M");
      double c  = cent();
      int index = profileIndex(multiplicityBins,c);
      
      //particle correlation block
      // Get trigger particles, charged hadrons  
      Particles trigg_h_Particles[PT_TRIGG_BINS];
      for (int ipt = 0; ipt < PT_TRIGG_BINS; ++ipt) {
        string pname = "APRIMTrigg" + toString(ipt);
        trigg_h_Particles[ipt] =
          applyProjection<ALICE::PrimaryParticles>(event,pname).particles();
      }
      // Get trigger particles, neutral hadrons  
      Particles trigg_V0_Particles[PT_TRIGG_BINS];
      for (int ipt = 0; ipt < PT_TRIGG_BINS; ++ipt) {
        string pname = "APRIMTrigg0" + toString(ipt);
        trigg_V0_Particles[ipt] =
          applyProjection<ALICE::PrimaryParticles>(event,pname).particles();
      }

      // Get associated particles particles, charged hadrons  
      Particles assocParticles[PT_ASSOC_BINS];
      for (int ipt = 0; ipt < PT_ASSOC_BINS; ++ipt) {
        string pname = "APRIMAssoc" + toString(ipt);
        assocParticles[ipt] =
          applyProjection<ALICE::PrimaryParticles>(event,pname).particles();
      }

      //trigger = any charged particle
      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg){
          if(trigg_h_Particles[ipt_trigg].size()==0)continue;
        for (const Particle& trigg : trigg_h_Particles[ipt_trigg]) {
          _counterChargedTriggers_mult[6][ipt_trigg]->fill(1);
          if(index>-1&&index<6)_counterChargedTriggers_mult[index][ipt_trigg]->fill(1);
          for (int ipt_assoc = 0; ipt_assoc < PT_ASSOC_BINS; ++ipt_assoc){
              if(assocParticles[ipt_assoc].size()==0)continue;
            for(const Particle& assoc : assocParticles[ipt_assoc]){
              if(assoc.pT() < trigg.pT()){
                fillbyparticles(_hist_hh_2D_ptassoc[ipt_trigg][ipt_assoc], trigg, assoc); 
                if(index>-1&&index<6)fillbyparticles(_hist_hh_2D_mult[index][ipt_trigg], trigg, assoc); 
                fillbyparticles(_hist_hh_2D_mult[6][ipt_trigg], trigg, assoc); 
              }
            }
          }
        }
      }

      //trigger = V0
      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg){
          if(trigg_V0_Particles[ipt_trigg].size()==0)continue;
        for (const Particle& triggV0 : trigg_V0_Particles[ipt_trigg]) {
          const int pid = abs(triggV0.pid());
          if(pid==310){
            _counterK0Triggers_mult[6][ipt_trigg]->fill(1);
            if(index>-1&&index<6)_counterK0Triggers_mult[index][ipt_trigg]->fill(1);
          }
          if(pid==3122){
            _counterLamTriggers_mult[6][ipt_trigg]->fill(1);
            if(index>-1&&index<6)_counterLamTriggers_mult[index][ipt_trigg]->fill(1);
          }
          for (int ipt_assoc = 0; ipt_assoc < PT_ASSOC_BINS; ++ipt_assoc){
              if(assocParticles[ipt_assoc].size()==0)continue;
            for(const Particle& assoc : assocParticles[ipt_assoc]){
              if(assoc.pT() < triggV0.pT() && pid==310){
                fillbyparticles(_hist_K0h_2D_mult[6][ipt_trigg], triggV0, assoc); 
                if(index>-1&&index<6)fillbyparticles(_hist_K0h_2D_mult[index][ipt_trigg], triggV0, assoc); 
                fillbyparticles(_hist_K0h_2D_ptassoc[ipt_trigg][ipt_assoc], triggV0, assoc); 
              }
              if(assoc.pT() < triggV0.pT() && pid==3122){
                fillbyparticles(_hist_Lamh_2D_mult[6][ipt_trigg], triggV0, assoc);  
                if(index>-1&&index<6)fillbyparticles(_hist_Lamh_2D_mult[index][ipt_trigg], triggV0, assoc);   
                fillbyparticles(_hist_Lamh_2D_ptassoc[ipt_trigg][ipt_assoc], triggV0, assoc);  
              }
            }
          }
        }
      }
      //end of particle correlation block
      
      //mixed event block
      const EventMixingFinalState& evmc = apply<EventMixingFinalState>(event, "EVMc");
      if (!evmc.hasMixingEvents()) return;

      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg){
          if(trigg_h_Particles[ipt_trigg].size()==0)continue;
        for (const Particle& trigg : trigg_h_Particles[ipt_trigg]) {
            if(evmc.particles().size()==0)continue;
          for (const Particle& assoc_mix : evmc.particles()){
            if (assoc_mix.pT() < trigg.pT())fillbyparticles(_hist_mix_hh, trigg, assoc_mix);   
          }
        }
      }
      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg){
          if(trigg_V0_Particles[ipt_trigg].size()==0)continue;
        for (const Particle& triggV0 : trigg_V0_Particles[ipt_trigg]) {
          const int pid = abs(triggV0.pid());
          if(evmc.particles().size()==0)continue;
          for (const Particle& assoc_mix : evmc.particles()){
            if(assoc_mix.pT() < triggV0.pT() && pid==310)fillbyparticles(_hist_mix_K0h, triggV0, assoc_mix); 
            if(assoc_mix.pT() < triggV0.pT() && pid==3122)fillbyparticles(_hist_mix_Lamh, triggV0, assoc_mix);   
          }
        }
      }
      //end of mixed event block
    }


    /// Finalize
    void finalize() {

      double mix_nomalisation_hh = (_hist_mix_hh.histo(0))->integral()/(_hist_mix_hh.histo(0))->numBins();
      double mix_scaling_hh[38];
      int i_mix=0;
      for(Histo1DPtr hist : _hist_mix_hh.histos()){
        mix_scaling_hh[i_mix]=(hist->integral()/hist->numBins())/mix_nomalisation_hh;
        i_mix++;
      }

      double mix_nomalisation_K0h = (_hist_mix_K0h.histo(0))->integral()/(_hist_mix_K0h.histo(0))->numBins();
      double mix_scaling_K0h[38];
      i_mix=0;
      for(Histo1DPtr hist : _hist_mix_K0h.histos()){
        mix_scaling_K0h[i_mix]=(hist->integral()/hist->numBins())/mix_nomalisation_K0h;
        i_mix++;
      }

      double mix_nomalisation_Lamh = (_hist_mix_Lamh.histo(0))->integral()/(_hist_mix_Lamh.histo(0))->numBins();
      double mix_scaling_Lamh[38];
      i_mix=0;
      for(Histo1DPtr hist : _hist_mix_Lamh.histos()){
        mix_scaling_Lamh[i_mix]=(hist->integral()/hist->numBins())/mix_nomalisation_Lamh;
        i_mix++;
      }
      
      for (int imult = 0; imult < MULT_BINS; ++imult) {
        for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg) {
          //cor. scaling + mixing
            if(_counterChargedTriggers_mult[imult][ipt_trigg]->sumW() >0)  for (Histo1DPtr hist : _hist_hh_2D_mult[imult][ipt_trigg].histos()) { scale(hist,1./_counterChargedTriggers_mult[imult][ipt_trigg]->sumW()); }
            i_mix=0;
            for (Histo1DPtr hist : _hist_hh_2D_mult[imult][ipt_trigg].histos()) { if(mix_scaling_hh[i_mix]>0) scale(hist,1./mix_scaling_hh[i_mix]); i_mix++; }
            if(_counterK0Triggers_mult[imult][ipt_trigg]->sumW() >0)  for (Histo1DPtr hist : _hist_K0h_2D_mult[imult][ipt_trigg].histos()) { scale(hist,1./_counterK0Triggers_mult[imult][ipt_trigg]->sumW()); }
            i_mix=0;
            for (Histo1DPtr hist : _hist_K0h_2D_mult[imult][ipt_trigg].histos()) { if(mix_scaling_K0h[i_mix]>0) scale(hist,1./mix_scaling_K0h[i_mix]); i_mix++; }
            if(_counterLamTriggers_mult[imult][ipt_trigg]->sumW() >0)  for (Histo1DPtr hist : _hist_Lamh_2D_mult[imult][ipt_trigg].histos()) { scale(hist,1./_counterLamTriggers_mult[imult][ipt_trigg]->sumW()); }
            i_mix=0;
            for (Histo1DPtr hist : _hist_Lamh_2D_mult[imult][ipt_trigg].histos()) { if(mix_scaling_Lamh[i_mix]>0) scale(hist,1./mix_scaling_Lamh[i_mix]); i_mix++; }
        
          //integration by eta
          S2DProjectionY(_hist_dPhi_hh_mult[imult][ipt_trigg], _hist_hh_2D_mult[imult][ipt_trigg]);
          S2DProjectionY(_hist_dPhi_K0h_mult[imult][ipt_trigg], _hist_K0h_2D_mult[imult][ipt_trigg]);
          S2DProjectionY(_hist_dPhi_Lamh_mult[imult][ipt_trigg], _hist_Lamh_2D_mult[imult][ipt_trigg]);
          //ZYAM
          ZYAM(_hist_dPhi_hh_mult_fin[imult][ipt_trigg], _hist_dPhi_hh_mult[imult][ipt_trigg]);
          ZYAM(_hist_dPhi_K0h_mult_fin[imult][ipt_trigg], _hist_dPhi_K0h_mult[imult][ipt_trigg]);
          ZYAM(_hist_dPhi_Lamh_mult_fin[imult][ipt_trigg], _hist_dPhi_Lamh_mult[imult][ipt_trigg]);
        }
      }

      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg) {
        for (int ipt_assoc = 0; ipt_assoc < PT_ASSOC_BINS; ++ipt_assoc) {
          //cor. scaling + mixing
          if(_counterChargedTriggers_mult[6][ipt_trigg]->sumW() >0)  for (Histo1DPtr hist : _hist_hh_2D_ptassoc[ipt_trigg][ipt_assoc].histos()) { scale(hist,1./_counterChargedTriggers_mult[6][ipt_trigg]->sumW()); }
          i_mix=0;
          for (Histo1DPtr hist : _hist_hh_2D_ptassoc[ipt_trigg][ipt_assoc].histos()) { if(mix_scaling_hh[i_mix]>0)scale(hist,1./mix_scaling_hh[i_mix]); i_mix++; }
          if(_counterK0Triggers_mult[6][ipt_trigg]->sumW() >0)  for (Histo1DPtr hist : _hist_K0h_2D_ptassoc[ipt_trigg][ipt_assoc].histos()) { scale(hist,1./_counterK0Triggers_mult[6][ipt_trigg]->sumW()); }
          i_mix=0;
          for (Histo1DPtr hist : _hist_K0h_2D_ptassoc[ipt_trigg][ipt_assoc].histos()) { if(mix_scaling_K0h[i_mix]>0)scale(hist,1./mix_scaling_K0h[i_mix]); i_mix++; }
          if(_counterLamTriggers_mult[6][ipt_trigg]->sumW() >0)  for (Histo1DPtr hist : _hist_Lamh_2D_ptassoc[ipt_trigg][ipt_assoc].histos()) { scale(hist,1./_counterLamTriggers_mult[6][ipt_trigg]->sumW()); }
          i_mix=0;
          for (Histo1DPtr hist : _hist_Lamh_2D_ptassoc[ipt_trigg][ipt_assoc].histos()) { if(mix_scaling_Lamh[i_mix]>0) scale(hist,1./mix_scaling_Lamh[i_mix]); i_mix++; }
        
          //integration by eta
          S2DProjectionY(_hist_dPhi_hh_ptassoc[ipt_trigg][ipt_assoc], _hist_hh_2D_ptassoc[ipt_trigg][ipt_assoc]);
          S2DProjectionY(_hist_dPhi_K0h_ptassoc[ipt_trigg][ipt_assoc], _hist_K0h_2D_ptassoc[ipt_trigg][ipt_assoc]);
          S2DProjectionY(_hist_dPhi_Lamh_ptassoc[ipt_trigg][ipt_assoc], _hist_Lamh_2D_ptassoc[ipt_trigg][ipt_assoc]);
          //Subtracting underlying event with the ZYAM method
          ZYAM(_hist_dPhi_hh_ptassoc_fin[ipt_trigg][ipt_assoc], _hist_dPhi_hh_ptassoc[ipt_trigg][ipt_assoc]);
          ZYAM(_hist_dPhi_K0h_ptassoc_fin[ipt_trigg][ipt_assoc], _hist_dPhi_K0h_ptassoc[ipt_trigg][ipt_assoc]);
          ZYAM(_hist_dPhi_Lamh_ptassoc_fin[ipt_trigg][ipt_assoc], _hist_dPhi_Lamh_ptassoc[ipt_trigg][ipt_assoc]);
        }
      }

      //PeakYield 
      pair<double,double> NearInterval={-0.9,0.9};
      pair<double,double> AwayInterval={M_PI-1.4,M_PI+1.4};

      for (int imult = MULT_BINS-1; imult > -1; --imult) { 
        //fig 3
        IntegratePeakByPT(_hist_hh_NearSideYield_mult[imult],_hist_dPhi_hh_mult_fin[imult], NearInterval,8);
        IntegratePeakByPT(_hist_K0h_NearSideYield_mult[imult],_hist_dPhi_K0h_mult_fin[imult], NearInterval,8);
        IntegratePeakByPT(_hist_Lamh_NearSideYield_mult[imult],_hist_dPhi_Lamh_mult_fin[imult], NearInterval,8);

        IntegratePeakByPT(_hist_hh_AwaySideYield_mult[imult],_hist_dPhi_hh_mult_fin[imult], AwayInterval,8);
        IntegratePeakByPT(_hist_K0h_AwaySideYield_mult[imult],_hist_dPhi_K0h_mult_fin[imult], AwayInterval,8);
        IntegratePeakByPT(_hist_Lamh_AwaySideYield_mult[imult],_hist_dPhi_Lamh_mult_fin[imult], AwayInterval,8);

        //fig 4
        if(imult<MULT_BINS-1){
          sdivide(_hist_hh_NearSideYield_mult[imult], _hist_hh_NearSideYield_mult[6], _hist_hh_NearSideYield_mult_ratio[imult],8,false);
          sdivide(_hist_hh_AwaySideYield_mult[imult], _hist_hh_AwaySideYield_mult[6], _hist_hh_AwaySideYield_mult_ratio[imult],8,false);

          sdivide(_hist_K0h_NearSideYield_mult[imult], _hist_K0h_NearSideYield_mult[6], _hist_K0h_NearSideYield_mult_ratio[imult],6,true);
          sdivide(_hist_K0h_AwaySideYield_mult[imult], _hist_K0h_AwaySideYield_mult[6], _hist_K0h_AwaySideYield_mult_ratio[imult],6,true);
          sdivide(_hist_Lamh_NearSideYield_mult[imult], _hist_Lamh_NearSideYield_mult[6], _hist_Lamh_NearSideYield_mult_ratio[imult],6,true);
          if(imult==1)sdivide(_hist_Lamh_AwaySideYield_mult[imult], _hist_Lamh_AwaySideYield_mult[6], _hist_Lamh_AwaySideYield_mult_ratio[imult],5,true);
          else sdivide(_hist_Lamh_AwaySideYield_mult[imult], _hist_Lamh_AwaySideYield_mult[6], _hist_Lamh_AwaySideYield_mult_ratio[imult],6,true);
        }
      }

      //fig 5
      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS; ++ipt_trigg)
      {
        IntegratePeakByPT(_hist_hh_NearSideYield_ptassoc[ipt_trigg], _hist_dPhi_hh_ptassoc_fin[ipt_trigg], NearInterval,ipt_trigg);
        IntegratePeakByPT(_hist_hh_AwaySideYield_ptassoc[ipt_trigg], _hist_dPhi_hh_ptassoc_fin[ipt_trigg], AwayInterval,ipt_trigg);
        IntegratePeakByPT(_hist_K0h_NearSideYield_ptassoc[ipt_trigg], _hist_dPhi_K0h_ptassoc_fin[ipt_trigg], NearInterval,ipt_trigg);
        IntegratePeakByPT(_hist_K0h_AwaySideYield_ptassoc[ipt_trigg], _hist_dPhi_K0h_ptassoc_fin[ipt_trigg], AwayInterval,ipt_trigg);
        if(ipt_trigg<7)IntegratePeakByPT(_hist_Lamh_NearSideYield_ptassoc[ipt_trigg], _hist_dPhi_Lamh_ptassoc_fin[ipt_trigg], NearInterval,ipt_trigg);
        
        if(ipt_trigg==6)IntegratePeakByPT(_hist_Lamh_AwaySideYield_ptassoc[ipt_trigg], _hist_dPhi_Lamh_ptassoc_fin[ipt_trigg], AwayInterval,2);
        else if(ipt_trigg==7) continue;
        else IntegratePeakByPT(_hist_Lamh_AwaySideYield_ptassoc[ipt_trigg], _hist_dPhi_Lamh_ptassoc_fin[ipt_trigg], AwayInterval,ipt_trigg);
      }

      //fig 8
      for (int imult = 0; imult < MULT_BINS; ++imult)
      {
        sdivide(_hist_K0h_NearSideYield_mult[imult], _hist_hh_NearSideYield_mult[imult], _hist_K0h_NearSideYield_ratioTo_hh_mult[imult],7,false);
        sdivide(_hist_K0h_AwaySideYield_mult[imult], _hist_hh_AwaySideYield_mult[imult], _hist_K0h_AwaySideYield_ratioTo_hh_mult[imult],7,false);
        
        sdivide(_hist_Lamh_NearSideYield_mult[imult], _hist_hh_NearSideYield_mult[imult], _hist_Lamh_NearSideYield_ratioTo_hh_mult[imult],7,false);
        sdivide(_hist_Lamh_AwaySideYield_mult[imult], _hist_hh_AwaySideYield_mult[imult], _hist_Lamh_AwaySideYield_ratioTo_hh_mult[imult],7,false);
      }

      //fig 9
      for (int ipt_trigg = 0; ipt_trigg < PT_TRIGG_BINS-1; ++ipt_trigg)
      {
        sdivide(_hist_K0h_NearSideYield_ptassoc[ipt_trigg], _hist_hh_NearSideYield_ptassoc[ipt_trigg], _hist_K0h_NearSideYield_ratioTo_hh_ptassoc[ipt_trigg],ipt_trigg,false);
        sdivide(_hist_K0h_AwaySideYield_ptassoc[ipt_trigg], _hist_hh_AwaySideYield_ptassoc[ipt_trigg], _hist_K0h_AwaySideYield_ratioTo_hh_ptassoc[ipt_trigg],ipt_trigg,false);
        
        sdivide(_hist_Lamh_NearSideYield_ptassoc[ipt_trigg], _hist_hh_NearSideYield_ptassoc[ipt_trigg], _hist_Lamh_NearSideYield_ratioTo_hh_ptassoc[ipt_trigg],ipt_trigg,false);
        if(ipt_trigg<6)sdivide(_hist_Lamh_AwaySideYield_ptassoc[ipt_trigg], _hist_hh_AwaySideYield_ptassoc[ipt_trigg], _hist_Lamh_AwaySideYield_ratioTo_hh_ptassoc[ipt_trigg],ipt_trigg,false);
        else sdivide(_hist_Lamh_AwaySideYield_ptassoc[ipt_trigg], _hist_hh_AwaySideYield_ptassoc[ipt_trigg], _hist_Lamh_AwaySideYield_ratioTo_hh_ptassoc[ipt_trigg],2,false);
      }
  }
    
    //@}


  private:

    static const int PT_TRIGG_BINS = 8;
    static const int PT_ASSOC_BINS = 10;
    static const int MULT_BINS = 7;
    static const int N_phibins = 72;
    /// @name Histograms
    //@{
    
    BinnedHistogram _hist_mix_hh, _hist_mix_K0h, _hist_mix_Lamh;

    CounterPtr _counterChargedTriggers_mult[MULT_BINS][PT_TRIGG_BINS];
    CounterPtr _counterK0Triggers_mult[MULT_BINS][PT_TRIGG_BINS];
    CounterPtr _counterLamTriggers_mult[MULT_BINS][PT_TRIGG_BINS];

    BinnedHistogram _hist_hh_2D_ptassoc[PT_TRIGG_BINS][PT_ASSOC_BINS];
    BinnedHistogram _hist_K0h_2D_ptassoc[PT_TRIGG_BINS][PT_ASSOC_BINS];
    BinnedHistogram _hist_Lamh_2D_ptassoc[PT_TRIGG_BINS][PT_ASSOC_BINS];    
    
    BinnedHistogram _hist_hh_2D_mult[MULT_BINS][PT_TRIGG_BINS];
    BinnedHistogram _hist_K0h_2D_mult[MULT_BINS][PT_TRIGG_BINS];
    BinnedHistogram _hist_Lamh_2D_mult[MULT_BINS][PT_TRIGG_BINS];
    
    Scatter2DPtr _hist_dPhi_hh_ptassoc[PT_TRIGG_BINS][PT_ASSOC_BINS];
    Scatter2DPtr _hist_dPhi_K0h_ptassoc[PT_TRIGG_BINS][PT_ASSOC_BINS];
    Scatter2DPtr _hist_dPhi_Lamh_ptassoc[PT_TRIGG_BINS][PT_ASSOC_BINS];
    
    Scatter2DPtr _hist_dPhi_hh_ptassoc_fin[PT_TRIGG_BINS][PT_ASSOC_BINS];
    Scatter2DPtr _hist_dPhi_K0h_ptassoc_fin[PT_TRIGG_BINS][PT_ASSOC_BINS];
    Scatter2DPtr _hist_dPhi_Lamh_ptassoc_fin[PT_TRIGG_BINS][PT_ASSOC_BINS];
    
    Scatter2DPtr _hist_dPhi_hh_mult[MULT_BINS][PT_TRIGG_BINS]; 
    Scatter2DPtr _hist_dPhi_K0h_mult[MULT_BINS][PT_TRIGG_BINS];
    Scatter2DPtr _hist_dPhi_Lamh_mult[MULT_BINS][PT_TRIGG_BINS];

    Scatter2DPtr _hist_dPhi_hh_mult_fin[MULT_BINS][PT_TRIGG_BINS];
    Scatter2DPtr _hist_dPhi_K0h_mult_fin[MULT_BINS][PT_TRIGG_BINS];
    Scatter2DPtr _hist_dPhi_Lamh_mult_fin[MULT_BINS][PT_TRIGG_BINS];

    Scatter2DPtr _hist_K0h_NearSideYield_ratioTo_hh_mult[MULT_BINS];
    Scatter2DPtr _hist_K0h_AwaySideYield_ratioTo_hh_mult[MULT_BINS];
    Scatter2DPtr _hist_Lamh_NearSideYield_ratioTo_hh_mult[MULT_BINS];
    Scatter2DPtr _hist_Lamh_AwaySideYield_ratioTo_hh_mult[MULT_BINS];

    Scatter2DPtr _hist_K0h_NearSideYield_ratioTo_hh_ptassoc[PT_TRIGG_BINS-1];
    Scatter2DPtr _hist_K0h_AwaySideYield_ratioTo_hh_ptassoc[PT_TRIGG_BINS-1];
    Scatter2DPtr _hist_Lamh_NearSideYield_ratioTo_hh_ptassoc[PT_TRIGG_BINS-1];
    Scatter2DPtr _hist_Lamh_AwaySideYield_ratioTo_hh_ptassoc[PT_TRIGG_BINS-1];

    Scatter2DPtr _hist_hh_NearSideYield_mult[MULT_BINS];
    Scatter2DPtr _hist_K0h_NearSideYield_mult[MULT_BINS];
    Scatter2DPtr _hist_Lamh_NearSideYield_mult[MULT_BINS];
    Scatter2DPtr _hist_hh_AwaySideYield_mult[MULT_BINS];
    Scatter2DPtr _hist_K0h_AwaySideYield_mult[MULT_BINS];
    Scatter2DPtr _hist_Lamh_AwaySideYield_mult[MULT_BINS];

    Scatter2DPtr _hist_hh_NearSideYield_mult_ratio[MULT_BINS];
    Scatter2DPtr _hist_K0h_NearSideYield_mult_ratio[MULT_BINS];
    Scatter2DPtr _hist_Lamh_NearSideYield_mult_ratio[MULT_BINS];
    Scatter2DPtr _hist_hh_AwaySideYield_mult_ratio[MULT_BINS];
    Scatter2DPtr _hist_K0h_AwaySideYield_mult_ratio[MULT_BINS];
    Scatter2DPtr _hist_Lamh_AwaySideYield_mult_ratio[MULT_BINS];

    Scatter2DPtr _hist_hh_NearSideYield_ptassoc[PT_TRIGG_BINS];
    Scatter2DPtr _hist_hh_AwaySideYield_ptassoc[PT_TRIGG_BINS];
    Scatter2DPtr _hist_K0h_NearSideYield_ptassoc[PT_TRIGG_BINS];
    Scatter2DPtr _hist_K0h_AwaySideYield_ptassoc[PT_TRIGG_BINS];
    Scatter2DPtr _hist_Lamh_NearSideYield_ptassoc[PT_TRIGG_BINS-1];
    Scatter2DPtr _hist_Lamh_AwaySideYield_ptassoc[PT_TRIGG_BINS-1];

    static const int Num_etabins = 39; 
    vector<double> multiplicityBins; 
    double etabins[Num_etabins];

    vector<double> bins_pt_trigg = {3. ,4. ,5. ,6. ,7. ,9. ,11. ,15. ,20.};
    vector<double> bins_pt_assoc = {1. ,2., 3. ,4. ,5. ,6. ,7. ,9. ,11. ,15., 20.};
    vector<string> mulsel_name = {"000-001", "001-003", "003-007", "007-015", "015-050", "050-100", "MB"};
    vector<string> bins_pt_trigg_name = {"3-4", "4-5", "5-6", "6-7", "7-9", "9-11", "11-15", "15-20"};
    vector<string> bins_pt_assoc_name = {"_1-2","_2-3", "_3-4", "_4-5", "_5-6", "_6-7", "_7-9", "_9-11", "_11-15", "_15-20"};
    
  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1891391);

}

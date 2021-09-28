// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  /// @brief Study of ordered hadron chains at 7 TeV
  class ATLAS_2017_I1624693 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1624693);

    /// @name Analysis methods
    //@{
    struct usedX {
      
      int locMin;
      int locMax;
      std::vector<std::pair<int,float> > chains; 
      
      // Constructor
      usedX(int min, int max, int ic, float mass) {
        locMin=min;
        locMax=max;
        chains.clear();
        chains.push_back(std::pair<int,float>(ic,mass));
      }
      
      // Constructor
      usedX(int min, int max) {
        locMin=min;
        locMax=max;
        chains.clear();
      }

      void add(int jc, float mass) {
        
        if (chains.size()) {
          std::vector<std::pair<int,float> >::iterator it=chains.begin();
          while ( it!=chains.end() && mass>(*it).second )  ++it;
          chains.insert(it,std::pair<int,float>(jc,mass));
        }
        else {
          chains.push_back(std::pair<int,float>(jc,mass));
        }
      }
    };


    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      ChargedFinalState cfs((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  0.1*GeV));
      declare(cfs,"CFS");

      // pion mass;
      pim = 0.1396;     

      /// @todo Book histograms here, e.g.:
      book(_DeltaQ , 1, 1, 1);
      book(_Delta3h, 2, 1, 1);
      book(_dalitz , 3, 1, 1);

      // auxiliary
      book(_h_nch, "_nch", 200, -0.5, 199.5);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //const double weight = event.weight();
      bool match =false;

      /// @todo Do the event by event analysis here
      const ChargedFinalState& had = applyProjection<ChargedFinalState>(event, "CFS");
      Particles hs=had.particles();
      int nch = hs.size(); 

      if (nch < 3)  return; 

      _h_nch->fill(1.*nch,1.);
     
      for (unsigned int i=0; i < hs.size() - 1; ++i) {
        for (unsigned int j=i+1; j < hs.size(); ++j) {
          double q12 = qq(hs[i],hs[j],match);
          if (match) _DeltaQ->fill(q12,-1.);
          else       _DeltaQ->fill(q12,1.);
        }
      }

      // chain selection
     
      std::vector<float> wchain;
      std::vector< std::vector<unsigned int> > rchains;
      std::vector< std::vector<float> > mchains;
      wchain.clear();  rchains.clear();  mchains.clear();
      for (unsigned int ip1 = 0; ip1< hs.size(); ++ip1 ) {
        wchain.push_back(1.);
        std::vector<unsigned int> cc(1,ip1);
        std::vector<float> mc;

        double qlmin=10000.; int ilmin=-1;
        for (unsigned ip2 = 0; ip2 < hs.size(); ++ip2) {
          if (ip2==ip1) continue;
          double ql = qq(hs[ip1],hs[ip2],match);
          if (!match) continue;    // looking for closest like-sign match      
          if (ql <qlmin) { qlmin=ql; ilmin=ip2;}
        }
        if (ilmin<0) { 
          wchain.back()=0.;
          mc.push_back(-1.);
        }
        else {     // search for unlike-sign match
          cc.push_back(ilmin);
          mc.push_back(qlmin);
          if (int(ip1)>ilmin && rchains[ilmin][1]==ip1) {
            // std::cout <<"exclusive match:"<< std::endl;
            wchain.back()=0.5; wchain[ilmin]=0.5;
          }
          
          double m3min=10000.; int ixmin=-1;
          for (unsigned ip2 = 0; ip2< hs.size(); ++ip2) {
            if (ip2==ip1 || int(ip2)==ilmin )  continue;
            double qx = qq(hs[ip1],hs[ip2],match);
            if (match)  continue;
            double qxl = qq(hs[ip2],hs[ilmin],match);
            double m3 = sqrt(9*pim*pim+qxl*qxl+qlmin*qlmin+qx*qx);
            if (m3 <m3min) { m3min=m3; ixmin=ip2;}
          }

          if (ixmin<0) {
            wchain.back()=0.;
            mc.push_back(-1.);
          }
          else {
            cc.push_back(ixmin);
            mc.push_back(m3min);
          }
        }
        rchains.push_back(cc);
        mchains.push_back(mc);
      }

      // cleanup: association rate for like-sign pairs should not exceed 2
      std::vector<float> assoc(hs.size(),0.);       // cache for association rate
      std::vector<bool> accept(rchains.size(), false);
      // loop over chains and accept lowest masses while watching the association rate
      int inext = 0;
      while ( inext>-1 ) {
        inext = -1; float cMin = 100000.;
        // find non-accepted chain with lowest Q_ls; dissolve chains if association count over 2 
        for (unsigned int ic=0; ic < rchains.size(); ++ic) {
          if (rchains[ic].size() < 2)  continue;
          if (accept[ic])  continue;
          if (mchains[ic][0] < cMin) { cMin = mchains[ic][0]; inext=ic; }
        }
        if (inext>-1 ) {
          unsigned int cloc0 = rchains[inext][0];
          unsigned int cloc1 = rchains[inext][1];
          if ( (assoc[cloc0] + 1. <= 2.)  &&  (assoc[cloc1] + 1. <= 2.) ) { // chain can be accepted
            accept[inext]=true;
            assoc[cloc0]+=1.;
            assoc[cloc1]+=1.;
            if (wchain[inext]==0.5) {  // accept the identical chain, too
              for (unsigned int ic=0; ic<hs.size(); ++ic) {
                if (rchains[ic][0] == cloc1 && rchains[ic][1] == cloc0) {   
                  accept[ic]=true;
                  break;
                }
              }
            }
          }
          else if ( assoc[cloc0]>1 ) { // association count filled up, discard chain
            accept[inext]=true;
            wchain[inext]=0.;
          }
          else {   // dissolve chain and find new association
            unsigned int i1 = rchains[inext][0];
            float mMn = 1000000.;
            int ipn = -1;
            for (unsigned int i2=0; i2<hs.size(); ++i2) {
              if (i1 == i2)  continue;
              double m = qq(hs[i1],hs[i2],match);
              if (!match)  continue;
              if (assoc[i2] > 1.)  continue;
              if (m > 0. && m <mMn ) { mMn = m; ipn = i2;}
            }
            if (ipn >= 0) {
              rchains[inext][1]=ipn;  mchains[inext][0]=mMn;
              // resolve chain weight : by default, it is 1.
              wchain[inext]=1.;
              // check exclusivity of pairing
              for (unsigned int ico=0; ico<hs.size(); ++ico) {
                if (int(rchains[ico][0]) == ipn && rchains[ico][1] == i1) {   // scale the contribution from both chains
                  wchain[ico]=0.5;
                  wchain[inext]=0.5;
                }
              }
              // add 3.member
              // continue with arbitrary match
              int ipnn=-1; float mMnn = 10000.;
              mMn = 1000000.;
              for (unsigned int ij=0; ij < hs.size(); ++ij) {
                rchains[inext].resize(2);
                float q02 = qq(hs[i1],hs[ij],match);
                if (match>0.)  continue;
                float q12 = qq(hs[ipn],hs[ij],match);
                double m3 = sqrt(9*pim*pim+q02*q02+mMn*mMn+q12*q12);
                if (m3>0. && m3 <mMnn ) { mMnn = m3; ipnn = ij; }
              }
              if (ipnn>=0) { rchains[inext].push_back(ipnn); rchains[inext][2]=ipnn; mchains[inext][1]=mMnn; }
              else {accept[inext]=true; wchain[inext]=0.;}
            }
            else { // chain not recovered
              wchain[inext]=0.;
              accept[inext]=true; 
            }
          }
        }
      }  // end loop over chains
 
      // cleanup: association rate for unlike-sign pairs
      // third member verification
      std::vector<bool> accept3(rchains.size(),false);
      // watch unlike-sign combinations used
      std::vector<usedX> used;
      // loop over chains and accept lowest masses while watching the association rate
      inext = 0;
      while ( inext>-1 ) {
        inext = -1; float cMin = 100000.;
        // find non-accepted chain with lowest mass; dissolve chains if association count over 3 
        for (unsigned int ic=0; ic < rchains.size(); ++ic) {
          if (rchains[ic].size() < 3 || !wchain[ic] || !accept[ic])  continue;
          if (accept3[ic])  continue;
          if (mchains[ic][1]<cMin) { cMin = mchains[ic][1]; inext=ic; }
        }
        // check association counts
        if (inext>-1 ) {
          unsigned int cloc0 = rchains[inext][0];
          unsigned int cloc1 = rchains[inext][1];
          unsigned int cloc2 = rchains[inext][2];

          // map use of unlike sign pairs
          int iu0 = -1; float w0=0.;
          for (unsigned int iu=0; iu<used.size(); ++iu) {
            if (fmin(cloc0,cloc2)==used[iu].locMin && fmax(cloc0,cloc2)==used[iu].locMax ) {
              iu0=iu;
              if (used[iu].chains.size() > 0) 
                for (unsigned int iw=0; iw<used[iu].chains.size(); ++iw)  w0+=wchain[used[iu].chains[iw].first];
              //used[iu].add(i1,mch[1]);
              break;
            }
          }
          if (iu0<0) { used.push_back(usedX(fmin(cloc0,cloc2),fmax(cloc0,cloc2)));iu0=used.size()-1; }
          int iu1 = -1; float w1=0.;
          for (unsigned int iu=0; iu<used.size(); ++iu) {
            if (fmin(cloc1,cloc2)==used[iu].locMin && fmax(cloc1,cloc2)==used[iu].locMax) {
              iu1=iu;
              if (used[iu].chains.size()>0) 
                for (unsigned int iw=0; iw<used[iu].chains.size(); iw++) w1 += wchain[used[iu].chains[iw].first];
              //used[iu].add(inext,mch[1]);
              break;
            }
          }
          if (iu1<0) { used.push_back(usedX(fmin(cloc1,cloc2),fmax(cloc1,cloc2))); iu1=used.size()-1; }

          if ( assoc[cloc2] < 3. && w0 < 2. && w1 < 2.) {
            accept3[inext] = true;
            assoc[cloc2] += 1.;
            used[iu0].add(inext, mchains[inext][1]);
            used[iu1].add(inext, mchains[inext][1]);
            if (wchain[inext]==0.5) {  // accept the identical chain, too
              for (unsigned int ic=0; ic< rchains.size(); ++ic) {
                if (rchains[ic][0]==cloc1 && rchains[ic][1] == cloc0) { 
                  accept3[ic]=true;
                  used[iu0].add(ic, mchains[ic][1]);
                  used[iu1].add(ic, mchains[ic][1]);
                  break;
                }
              }
            }
          }
          else { // find new association
            int i1 = rchains[inext][0];
            int i2 = rchains[inext][1];
            float mMn = 1000000.;
            int ipn=-1; int iploc=-1;
            rchains[inext].pop_back();
            for (unsigned int i3 = 0; i3 < hs.size(); ++i3) {
              double q02 = qq(hs[i1],hs[i3],match);
              if (match > 0.)  continue;
              if (assoc[i3] > 3-wchain[inext])  continue;
              // check pair association
              w0=0.; w1=0.;
              for (unsigned int iu=0; iu<used.size(); ++iu) {
                if (fmin(cloc0,i3)==used[iu].locMin && fmax(cloc0,i3)==used[iu].locMax ) {
                  if (used[iu].chains.size() > 0)
                    for (unsigned int iw=0; iw<used[iu].chains.size(); ++iw)  w0 += wchain[used[iu].chains[iw].first];
                }
                if (fmin(cloc1,i3)==used[iu].locMin && fmax(cloc1,i3)==used[iu].locMax ) {
                  if (used[iu].chains.size()>0) 
                    for (unsigned int iw=0; iw<used[iu].chains.size(); ++iw)  w1 += wchain[used[iu].chains[iw].first];
                }
              }
              if (w0+wchain[inext]>2. || w1+wchain[inext]>2.) continue;

              float q12 = qq(hs[i2],hs[i3],match);
              float q01 = qq(hs[i1],hs[i2],match);
              float m = sqrt(9*pim*pim+q02*q02+q01*q01+q12*q12);
              if (m>0. && m <mMn ) { mMn = m; ipn = i3; iploc = i3; }
            }
            if (ipn>=0) { 
              rchains[inext].push_back(ipn); rchains[inext][2]=iploc; mchains[inext][1]=mMn;
            }
            else { // chain not recovered
              wchain[inext]=0.;
            }
          }
        }
      }  // end loop over chains
      // end 3rd member optimization

      for (unsigned int ip=0; ip < wchain.size(); ++ip) {
        if (!wchain[ip])  continue;
        if (rchains[ip].size() < 3)  continue;
        float m3min = mchains[ip][1];
        if (m3min > 0.59)  continue;
        // dalitz plot
        std::pair<float,float> dd = dalitz3(hs[rchains[ip][0]], hs[rchains[ip][1]], hs[rchains[ip][2]]);
        _dalitz->fill(dd.first,dd.second,1.*wchain[ip]);
        // Delta(Q) spectra
        float qlmin = mchains[ip][0]; 
        float qxmin = qq(hs[rchains[ip][0]], hs[rchains[ip][2]], match); 
        float xlmin = qq(hs[rchains[ip][1]], hs[rchains[ip][2]], match); 
        _Delta3h->fill(qxmin, 0.5*wchain[ip]);
        _Delta3h->fill(xlmin, 0.5*wchain[ip]);
        _Delta3h->fill(qlmin, -1.*wchain[ip]);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      
      // normalize by the number of charged particles
      // counter automatic division by bin size 
      double norm = 0.01 / (_h_nch->xMean()*_h_nch->numEntries());
      _dalitz->scaleW(norm);
      _DeltaQ->scaleW(norm);
      _Delta3h->scaleW(norm);
      
    }

    //@}
    double qq(const Particle& gp1, const Particle& gp2, bool& match) {
      match = gp1.charge() * gp2.charge() > 0;
      FourMomentum p1, p2;
      p1.setPM(gp1.px(), gp1.py(), gp1.pz(), pim);
      p2.setPM(gp2.px(), gp2.py(), gp2.pz(), pim);
      return sqrt(fmax(0., (p1 + p2).mass2() - 4*pim*pim));
    }

    std::pair<float,float> dalitz3(const Particle& gp1, const Particle& gp2, const Particle& gp3) const {
      
      float p1= gp1.pt();
      float p2= gp2.pt();
      float p3= gp3.pt();
      float th1 = gp1.theta();
      float th2 = gp2.theta();
      float th3 = gp3.theta();
      float ph1 = gp1.phi();
      float ph2 = gp2.phi();
      float ph3 = gp3.phi();
      float e1 = sqrt(p1*p1+pim*pim);
      float e2 = sqrt(p2*p2+pim*pim);
      float e3 = sqrt(p3*p3+pim*pim);
      
      float p1x = p1*cos(ph1)*sin(th1);
      float p1y = p1*sin(ph1)*sin(th1);
      float p1z = p1*cos(th1);
      
      float p2x = p2*cos(ph2)*sin(th2);
      float p2y = p2*sin(ph2)*sin(th2);
      float p2z = p2*cos(th2);
      
      float p3x = p3*cos(ph3)*sin(th3);
      float p3y = p3*sin(ph3)*sin(th3);
      float p3z = p3*cos(th3);
      
      float px = p1x+p2x+p3x;
      float py = p1y+p2y+p3y;
      float pz = p1z+p2z+p3z;
      float ap = sqrt(px*px+py*py+pz*pz);
      float e=e1+e2+e3; 
      
      float beta = ap/e;
      float gamma = 1./sqrt(1-beta*beta);
      
      float p1l = (p1x*px+p1y*py+p1z*pz)/ap;
      float p2l = (p2x*px+p2y*py+p2z*pz)/ap;
      float p3l = (p3x*px+p3y*py+p3z*pz)/ap;
      
      float e1_boost = gamma*e1-gamma*beta*p1l;
      float e2_boost = gamma*e2-gamma*beta*p2l;
      float e3_boost = gamma*e3-gamma*beta*p3l;
      
      float Q = sqrt(e*e-ap*ap)-3*pim;
      
      return std::pair<float,float>(sqrt(3.)*(e1_boost-e2_boost)/Q , 3*(e3_boost-pim)/Q-1.);
    }
    
  private:

    // Data members like post-cuts event weight counters go here
    float pim;

  private:

    /// @name Histograms
    Histo1DPtr  _DeltaQ;
    Histo1DPtr _Delta3h;
    Histo1DPtr   _h_nch;
    Histo2DPtr _dalitz;
    //@}
  };

  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1624693);
}

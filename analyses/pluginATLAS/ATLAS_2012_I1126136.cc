// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2012_I1126136 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1126136()
      : Analysis("ATLAS_2012_I1126136")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialize projections before the run
    void init() {

      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 20*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");

      // Jet finder
      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.9),"vfs");

      // Book histograms
      book(_count_SR_A     ,"count_SR_A"    , 1, 0., 1.);
      book(_count_SR_B     ,"count_SR_B"    , 1, 0., 1.);

      book(_hist_mjjj1  ,"hist_mjjj1" , 30 , 0.   , 600.  );
      book(_hist_mjjj2  ,"hist_mjjj2" , 30 , 0.   , 600.  );
      book(_hist_ETmiss ,"hist_ETmiss", 20 , 100. , 600.  );
      book(_hist_mT2    ,"hist_mT2"   , 200, 0.   , 1000. );

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // pTmiss
      FourMomentum pTmiss;
      for (const Particle& p : apply<VisibleFinalState>(event, "vfs").particles() ) {
        pTmiss -= p.momentum();
      }
      double ETmiss = pTmiss.pT();

      // require eTmiss > 150
      if (ETmiss < 150*GeV) vetoEvent;

      // get the candiate jets
      Jets cand_jets;
      for ( const Jet& jet : apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if (jet.abseta() < 4.5) cand_jets.push_back(jet);
      }

      // find the electrons
      Particles cand_e;
      for( const Particle& e : apply<IdentifiedFinalState>(event, "elecs").particlesByPt()) {
        // remove any leptons within 0.4 of any candidate jets
        bool e_near_jet = false;
        for ( const Jet& jet : cand_jets ) {
          double dR = deltaR(e, jet);
          if (inRange(dR, 0.2, 0.4)) {
            e_near_jet = true;
            break;
          }
        }
        if ( e_near_jet ) continue;
        cand_e.push_back(e);
      }

      // find the muons
      Particles cand_mu;
      for( const Particle& mu : apply<IdentifiedFinalState>(event, "muons").particlesByPt()) {
        // remove any leptons within 0.4 of any candidate jets
        bool mu_near_jet = false;
        for ( const Jet& jet : cand_jets ) {
          if ( deltaR(mu, jet) < 0.4 ) {
            mu_near_jet = true;
            break;
          }
        }
        if ( mu_near_jet ) continue;
        cand_mu.push_back(mu);
      }

      // veto events with leptons
      if( ! cand_e.empty() || ! cand_mu.empty() )
        vetoEvent;

      // discard jets that overlap with electrons
      Jets recon_jets;
      for ( const Jet& jet : cand_jets ) {
        if (jet.abseta() > 2.8 || jet.pT() < 30*GeV) continue;
        bool away_from_e = true;
        for (const Particle& e : cand_e ) {
          if ( deltaR(e, jet) < 0.2 ) {
            away_from_e = false;
            break;
          }
        }
        if ( away_from_e ) recon_jets.push_back( jet );
      }

      // find b jets
      Jets tight_bjets,loose_bjets;
      for(const Jet& jet : recon_jets) {
        /// @todo Should be abseta?
        if (!jet.bTagged() && jet.eta()>2.5) continue;
        double prob = rand()/static_cast<double>(RAND_MAX);
        if (prob <= 0.60) tight_bjets.push_back(jet);
        if (prob <= 0.75) loose_bjets.push_back(jet);
      }

      // require >=1 tight or >=2 loose b-jets
      if (! ( !tight_bjets.empty() || loose_bjets.size()>=2) )
        vetoEvent;

      // must be at least 6 jets with pT>30
      if (recon_jets.size()<6 ) vetoEvent;

      // hardest > 130
      if (recon_jets[0].perp() < 130. ) vetoEvent;

      // three hardest jets must be separated from etmiss
      for (unsigned int ix=0;ix<3;++ix) {
        if (deltaPhi(recon_jets[ix].momentum(),pTmiss)<0.2*PI)
          vetoEvent;
      }

      // remove events with tau like jets
      for (unsigned int ix=3;ix<recon_jets.size();++ix) {
        // skip jets seperated from eTmiss
        if (deltaPhi(recon_jets[ix].momentum(),pTmiss)>=0.2*PI)
          continue;
        // check the number of tracks between 1 and 4
        unsigned int ncharged=0;
        for ( const Particle & particle : recon_jets[ix].particles()) {
          if (PID::charge3(particle.pid())!=0) ++ncharged;
        }
        if (ncharged==0 || ncharged>4) continue;
        // calculate transverse mass and reject if < 100
        double mT = 2.*recon_jets[ix].perp()*ETmiss
          -recon_jets[ix].px()*pTmiss.px()
          -recon_jets[ix].py()*pTmiss.py();
        if (mT<100.) vetoEvent;
      }

      // if 2 loose b-jets apply mT cut
      if (loose_bjets.size()>=2) {
        // find b-jet closest to eTmiss
        double minR(1e30);
        unsigned int ijet(0);
        for(unsigned int ix=0;ix<loose_bjets.size();++ix) {
          double dR = deltaR(loose_bjets[ix].momentum(),pTmiss);
          if(dR<minR) {
            minR=dR;
            ijet = ix;
          }
        }
        double  mT = 2.*loose_bjets[ijet].perp()*ETmiss
          -loose_bjets[ijet].px()*pTmiss.px()
          -loose_bjets[ijet].py()*pTmiss.py();
        if(mT<170.) vetoEvent;
      }

      // 1 tight b-jet apply mT cut
      if(tight_bjets.size()==1) {
        for(unsigned int ix=0;ix<4;++ix) {
          double mT = 2.*recon_jets[ix].perp()*ETmiss
            -recon_jets[ix].px()*pTmiss.px()
            -recon_jets[ix].py()*pTmiss.py();
          if(mT<175.) vetoEvent;
        }
      }

      // find the closest triplet of jets in (eta,phi)
      unsigned int j1(0),j2(0),j3(0);
      double minR2(1e30);
      for(unsigned int i1=0;i1<recon_jets.size();++i1) {
        for(unsigned int i2=i1+1;i2<recon_jets.size();++i2) {
          for(unsigned int i3=i2+1;i3<recon_jets.size();++i3) {
            double delR2 =
              sqr(deltaR(recon_jets[i1].momentum(),recon_jets[i2].momentum())) +
              sqr(deltaR(recon_jets[i1].momentum(),recon_jets[i3].momentum())) +
              sqr(deltaR(recon_jets[i2].momentum(),recon_jets[i3].momentum()));
            if(delR2<minR2) {
              minR2=delR2;
              j1=i1;
              j2=i2;
              j3=i3;
            }
          }
        }
      }
      // 4-momentum and mass of first triplet
      FourMomentum pjjj1 = recon_jets[j1].momentum() +
        recon_jets[j2].momentum()+ recon_jets[j3].momentum();
      double mjjj1 = pjjj1.mass();

      // find the second triplet
      unsigned int j4(0),j5(0),j6(0);
      minR2=0.;
      for(unsigned int i1=0;i1<recon_jets.size();++i1) {
        if(i1==j1||i1==j2||i1==j3) continue;
        for(unsigned int i2=i1+1;i2<recon_jets.size();++i2) {
          if(i2==j1||i2==j2||i2==j3) continue;
          for(unsigned int i3=i2+1;i3<recon_jets.size();++i3) {
            if(i3==j1||i3==j2||i3==j3) continue;
            double delR2 =
              sqr(deltaR(recon_jets[i1].momentum(),recon_jets[i2].momentum())) +
              sqr(deltaR(recon_jets[i1].momentum(),recon_jets[i3].momentum())) +
              sqr(deltaR(recon_jets[i2].momentum(),recon_jets[i3].momentum()));
            if(delR2<minR2) {
              minR2=delR2;
              j4=i1;
              j5=i2;
              j6=i3;
            }
          }
        }
      }

      // 4-momentum and mass of first triplet
      FourMomentum pjjj2 = recon_jets[j4].momentum() +
        recon_jets[j5].momentum()+ recon_jets[j6].momentum();
      double mjjj2 = pjjj2.mass();

      _hist_mjjj1->fill(mjjj1,weight);
      _hist_mjjj2->fill(mjjj2,weight);
      // require triplets in 80<mjjj<270
      if(mjjj1<80.||mjjj1>270.||mjjj2<80.||mjjj2>270.)
        vetoEvent;

      // counts in signal regions
      _count_SR_A->fill(0.5,weight);
      if(ETmiss>260.) _count_SR_B->fill(0.5,weight);

      _hist_ETmiss->fill(ETmiss,weight);
      const double m_T2 = mT2(pjjj1, pjjj2, pTmiss, 0.0); // zero mass invisibles
      _hist_mT2->fill(m_T2,weight);
    }
    //@}


    void finalize() {

      double norm = 4.7* crossSection()/sumOfWeights()/femtobarn;
      scale(_count_SR_A ,     norm );
      scale(_count_SR_B ,     norm );
      scale(_hist_mjjj1 , 20.*norm );
      scale(_hist_ETmiss, 50.*norm );
      scale(_hist_mjjj2 , 20.*norm );
      scale(_hist_mT2   ,     norm );

    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _count_SR_A;
    Histo1DPtr _count_SR_B;

    Histo1DPtr _hist_mjjj1;
    Histo1DPtr _hist_mjjj2;
    Histo1DPtr _hist_ETmiss;
    Histo1DPtr _hist_mT2;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1126136);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

namespace H1_2007_I746380_PROJECTIONS {
  
  /// Projection to find the largest gaps and the masses of the two
  /// systems separated by the gap. Based on the HZTools gap-finding
  /// method (hzhadgap.F). Note that gaps are found in the HCM frame.
  ///
  /// @author Christine O. Rasmussen.
  class RapidityGap : public Projection {

  public:

    /// Type of DIS boost to apply
    enum Frame { HCM, LAB, XCM };

    RapidityGap() {
      setName("RapidityGap");
      declare(DISKinematics(), "DISKIN");
      declare(DISFinalState(DISFinalState::BoostFrame::HCM), "DISFS");
    }

    DEFAULT_RIVET_PROJ_CLONE(RapidityGap);

    const double M2X()               const {return _M2X;}
    const double M2Y()               const {return _M2Y;}
    const double t()                 const {return _t;}
    const double gap()               const {return _gap;}
    const double gapUpp()            const {return _gapUpp;}
    const double gapLow()            const {return _gapLow;}
    const double EpPzX(Frame f) const {
      if (f == LAB) return _ePpzX_LAB;
      else if (f == XCM) return _ePpzX_XCM;
      else return _ePpzX_HCM;
    }
    const double EmPzX(Frame f) const {
      if (f == LAB) return _eMpzX_LAB;
      else if (f == XCM) return _eMpzX_XCM;
      else return _eMpzX_HCM;
    }
    const FourMomentum pX(Frame f) const {
      if (f == LAB) return _momX_LAB;
      else if (f == XCM) return _momX_XCM;
      else return _momX_HCM;
    }
    const FourMomentum pY(Frame f) const {
      if (f == LAB) return _momY_LAB;
      else if (f == XCM) return _momY_XCM;
      else return _momY_HCM;
    }
    const Particles& systemX(Frame f) const {
      if (f == LAB) return _pX_LAB;
      else if (f == XCM) return _pX_XCM;
      else return _pX_HCM;
    }
    const Particles& systemY(Frame f) const {
      if (f == LAB) return _pY_LAB;
      else if (f == XCM) return _pY_XCM;
      else return _pY_HCM;
    }

  protected:

    virtual CmpState compare(const Projection& p) const {
      const RapidityGap& other = pcast<RapidityGap>(p);
      return mkNamedPCmp(other, "DISKIN") || mkNamedPCmp(other, "DISFS");
    }

    virtual void project(const Event& e){
      const DISKinematics& dk = apply<DISKinematics>(e, "DISKIN");
      const Particles& p      = apply<DISFinalState>(e, "DISFS").particles(cmpMomByEta);
      findgap(p, dk);
    }

    void clearAll(){
      _M2X = _M2Y = _t = _gap = 0.;
      _gapUpp = _gapLow = -8.;
      _ePpzX_HCM = _eMpzX_HCM =_ePpzX_LAB = _eMpzX_LAB = _ePpzX_XCM = _eMpzX_XCM = 0.;
      _momX_HCM.setPE(0., 0., 0., 0.);
      _momY_HCM.setPE(0., 0., 0., 0.);
      _momX_XCM.setPE(0., 0., 0., 0.);
      _momY_XCM.setPE(0., 0., 0., 0.);
      _momX_LAB.setPE(0., 0., 0., 0.);
      _momY_LAB.setPE(0., 0., 0., 0.);
      _pX_HCM.clear();
      _pY_HCM.clear();
      _pX_XCM.clear();
      _pY_XCM.clear();
      _pX_LAB.clear();
      _pY_LAB.clear();
    }

    void findgap(const Particles& particles, const DISKinematics& diskin){

      clearAll();

      // Begin by finding largest gap and gapedges between all final
      // state particles in HCM frame.
      int nP  = particles.size();
      int dir = diskin.orientation();
      for (int i = 0; i < nP-1; ++i){
        double tmpGap = abs(particles[i+1].eta() - particles[i].eta());
        if (tmpGap > _gap) {
          _gap    = tmpGap;
          _gapLow = (dir > 0) ? particles[i].eta() : dir * particles[i+1].eta();
          _gapUpp = (dir > 0) ? particles[i+1].eta() : dir * particles[i].eta();
        }
      }

      // Define the two systems X and Y.
      Particles tmp_pX, tmp_pY;
      for (const Particle& ip : particles) {
        if (dir * ip.eta() > _gapLow) tmp_pX.push_back(ip);
        else tmp_pY.push_back(ip);
      }

      Particles pX, pY;
      pX = (dir < 0) ? tmp_pY : tmp_pX;
      pY = (dir < 0) ? tmp_pX : tmp_pY;

      // Find variables related to HCM frame.
      // Note that HCM has photon along +z, as opposed to
      // H1 where proton is along +z. This results in a sign change
      // as compared to H1 papers!

      // X - side
      FourMomentum momX;
      for (const Particle& jp : pX) {
        momX  += jp.momentum();
        _ePpzX_HCM += jp.E() - jp.pz(); // Sign + => -
        _eMpzX_HCM += jp.E() + jp.pz(); // Sign - => +
      }
      _momX_HCM = momX;
      _pX_HCM   = pX;
      _M2X      = _momX_HCM.mass2();

      // Y - side
      FourMomentum momY;
      for (const Particle& kp : pY) momY += kp.momentum();
      _momY_HCM = momY;
      _pY_HCM   = pY;
      _M2Y      = _momY_HCM.mass2();

      // Find variables related to LAB frame
      const LorentzTransform hcmboost   = diskin.boostHCM();
      const LorentzTransform hcminverse = hcmboost.inverse();
      _momX_LAB = hcminverse.transform(_momX_HCM);
      _momY_LAB = hcminverse.transform(_momY_HCM);

      // Find momenta in XCM frame. Note that it is HCM frame that is
      // boosted, resulting in a sign change later!
      const bool doXCM = (momX.betaVec().mod2() < 1.);
      if (doXCM) {
        const LorentzTransform xcmboost =
          LorentzTransform::mkFrameTransformFromBeta(momX.betaVec());
        _momX_XCM = xcmboost.transform(momX);
        _momY_XCM = xcmboost.transform(momY);
      }

      for (const Particle& jp : pX) {
        // Boost from HCM to LAB.
        FourMomentum lab = hcminverse.transform(jp.momentum());
        _ePpzX_LAB += lab.E() + dir * lab.pz();
        _eMpzX_LAB += lab.E() - dir * lab.pz();
        Particle plab = jp;
        plab.setMomentum(lab);
        _pX_LAB.push_back(plab);
        // Set XCM. Note that since HCM frame is boosted to XCM frame,
        // we have a sign change
        if (doXCM) {
          const LorentzTransform xcmboost =
            LorentzTransform::mkFrameTransformFromBeta(_momX_HCM.betaVec());
          FourMomentum xcm = xcmboost.transform(jp.momentum());
          _ePpzX_XCM += xcm.E() - xcm.pz(); // Sign + => -
          _eMpzX_XCM += xcm.E() + xcm.pz(); // Sign - => +
          Particle pxcm = jp;
          pxcm.setMomentum(xcm);
          _pX_XCM.push_back(pxcm);
        }
      }

      for (const Particle& jp : pY) {
        // Boost from HCM to LAB
        FourMomentum lab = hcminverse.transform(jp.momentum());
        Particle plab = jp;
        plab.setMomentum(lab);
        _pY_LAB.push_back(plab);
        // Boost from HCM to XCM
        if (doXCM) {
          const LorentzTransform xcmboost =
            LorentzTransform::mkFrameTransformFromBeta(_momX_HCM.betaVec());
          FourMomentum xcm = xcmboost.transform(jp.momentum());
          Particle pxcm = jp;
          pxcm.setMomentum(xcm);
          _pY_XCM.push_back(pxcm);
        }
      }

      // Find t: Currently can only handle gap on proton side.
      // @TODO: Expand to also handle gap on photon side
      // Boost p from LAB to HCM frame to find t.
      const FourMomentum proton = hcmboost.transform(diskin.beamHadron().momentum());
      FourMomentum pPom         = proton - _momY_HCM;
      _t                        = pPom * pPom;

    }

  private:

    double _M2X, _M2Y, _t;
    double _gap, _gapUpp, _gapLow;
    double _ePpzX_LAB, _eMpzX_LAB, _ePpzX_HCM, _eMpzX_HCM, _ePpzX_XCM, _eMpzX_XCM;
    FourMomentum _momX_HCM, _momY_HCM,_momX_LAB, _momY_LAB, _momX_XCM, _momY_XCM;
    Particles _pX_HCM, _pY_HCM, _pX_LAB, _pY_LAB, _pX_XCM, _pY_XCM;

  };

  /// Projection to boost system X (photon+Pomeron) particles into its rest frame.
  ///
  /// @author Ilkka Helenius
  class BoostedXSystem : public FinalState {
  public:

    BoostedXSystem(const FinalState& fs) {
      setName("BoostedXSystem");
      declare(fs,"FS");
      declare(RapidityGap(), "RAPGAP");
    }

    // Return the boost to XCM frame.
    const LorentzTransform& boost() const { return _boost; }

    DEFAULT_RIVET_PROJ_CLONE(BoostedXSystem);

  protected:

    // Apply the projection on the supplied event.
    void project(const Event& e){

      const RapidityGap& rg = apply<RapidityGap>(e, "RAPGAP");

      // Total momentum of the system X.
      const FourMomentum pX = rg.pX(RapidityGap::HCM);

      // Reset the boost. Is there a separate method for this?
      _boost = combine(_boost, _boost.inverse());

      // Define boost only when numerically safe, otherwise negligible.
      if (pX.betaVec().mod2() < 1.)
        _boost = LorentzTransform::mkFrameTransformFromBeta(pX.betaVec());

      // Boost the particles from system X.
      _theParticles.clear();
      _theParticles.reserve(rg.systemX(RapidityGap::HCM).size());
      for (const Particle& p : rg.systemX(RapidityGap::HCM)) {
        Particle temp = p;
        temp.setMomentum(_boost.transform(temp.momentum()));
        _theParticles.push_back(temp);
      }

    }

    // Compare projections.
    CmpState compare(const Projection& p) const {
      const BoostedXSystem& other = pcast<BoostedXSystem>(p);
      return mkNamedPCmp(other, "RAPGAP") || mkNamedPCmp(other, "FS");
    }

  private:

    LorentzTransform _boost;

  };

  }



  /// @brief H1 diffractive dijets
  ///
  /// Diffractive dijets H1 with 920 GeV p and 27.5 GeV e
  /// Note tagged protons!
  ///
  /// @author Christine O. Rasmussen
  class H1_2007_I746380 : public Analysis {
  public:

    typedef H1_2007_I746380_PROJECTIONS::RapidityGap RapidityGap;
    typedef H1_2007_I746380_PROJECTIONS::BoostedXSystem BoostedXSystem;

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(H1_2007_I746380);

    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {

      declare(DISKinematics(), "Kinematics");
      const DISFinalState& disfs = declare(DISFinalState(DISFinalState::BoostFrame::HCM), "DISFS");
      const BoostedXSystem& disfsXcm = declare( BoostedXSystem(disfs), "BoostedXFS");
      declare(FastJets(disfsXcm, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 1.0,
                       JetAlg::Muons::ALL, JetAlg::Invisibles::NONE, nullptr), "DISFSJets");
      declare(RapidityGap(), "RapidityGap");

      // Book histograms from REF data
      book(_h_DIS_dsigdzPom, 1, 1, 1);
      book(_h_DIS_dsigdlogXpom, 2, 1, 1);
      book(_h_DIS_dsigdW, 3, 1, 1);
      book(_h_DIS_dsigdQ2, 4, 1, 1);
      book(_h_DIS_dsigdEtJet1, 5, 1, 1);
      book(_h_DIS_dsigdAvgEta, 6, 1, 1);
      book(_h_DIS_dsigdDeltaEta, 7, 1, 1);

      book(_h_PHO_dsigdzPom, 8, 1, 1);
      book(_h_PHO_dsigdxGam, 9, 1, 1);
      book(_h_PHO_dsigdlogXpom, 10, 1, 1);
      book(_h_PHO_dsigdW, 11, 1, 1);
      book(_h_PHO_dsigdEtJet1, 12, 1, 1);
      book(_h_PHO_dsigdAvgEta, 13, 1, 1);
      book(_h_PHO_dsigdDeltaEta, 14, 1, 1);
      book(_h_PHO_dsigdMjets, 15, 1, 1);

      isDIS  = false;
      nVeto0 = 0;
      nVeto1 = 0;
      nVeto2 = 0;
      nVeto3 = 0;
      nVeto4 = 0;
      nVeto5 = 0;
      nPHO   = 0;
      nDIS   = 0;
    }

    // Do the analysis
    void analyze(const Event& event) {

      // Event weight
      isDIS  = false;

      // Projections - special handling of events where no proton found:
      const RapidityGap&    rg = apply<RapidityGap>(event, "RapidityGap");
      const DISKinematics& kin = apply<DISKinematics>(event, "Kinematics");
      const BoostedXSystem& disfsXcm = apply<BoostedXSystem>( event, "BoostedXFS");

      // Determine kinematics: H1 has +z = proton direction
      int dir   = kin.orientation();
      double W2 = kin.W2();
      double W  = sqrt(W2);
      double y  = kin.y();
      double Q2 = kin.Q2();

      // Separate into DIS and PHO regimes else veto
      if (!inRange(W, 165.*GeV, 242.*GeV)) vetoEvent;
      if (Q2 < 0.01*GeV2) {
        isDIS = false;
        ++nPHO;
      } else if (inRange(Q2, 4.0*GeV2, 80.*GeV2)) {
        isDIS = true;
        ++nDIS;
      } else {
        vetoEvent;
      }
      ++nVeto0;

      // Find diffractive variables as defined in paper.
      const double M2Y  = rg.M2Y();
      const double M2X  = rg.M2X();
      const double abst = abs(rg.t());
      const double xPom = (isDIS) ? (Q2 + M2X) / (Q2 + W2) :
                          rg.EpPzX(RapidityGap::LAB) / (2. * kin.beamHadron().E());

      // Veto if outside allowed region
      if (sqrt(M2Y) > 1.6*GeV)    vetoEvent;
      ++nVeto1;
      if (abst > 1.0*GeV2) vetoEvent;
      ++nVeto2;
      if (xPom > 0.03)     vetoEvent;
      ++nVeto3;

      // Jet selection. Note jets are found in photon-proton (XCM)
      // frame, but eta cut is applied in lab frame!
      Cut jetcuts = Cuts::Et > 4.* GeV;
      Jets jets   = apply<FastJets>(event, "DISFSJets").jets(jetcuts, cmpMomByEt);
      // Veto if not dijets and if Et_j1 < 5.0
      if (jets.size() < 2)       vetoEvent;
      if (jets[0].Et() < 5.*GeV) vetoEvent;
      ++nVeto4;
      // Find Et_jet1 and deltaEta* in XCM frame
      double EtJet1       = jets[0].Et() * GeV;
      double etaXCMJet1   = jets[0].eta();
      double etaXCMJet2   = jets[1].eta();
      double deltaEtaJets = abs(etaXCMJet1 - etaXCMJet2);

      // Transform from XCM to HCM
      const LorentzTransform xcmboost = disfsXcm.boost();
      for (int i = 0; i < 2; ++i) jets[i].transformBy(xcmboost.inverse());
      // Find mass of jets and EpPz, EmPz of jets
      FourMomentum momJets = jets[0].momentum() + jets[1].momentum();
      double M2jets   = momJets.mass2();
      double EpPzJets = 0.;
      double EmPzJets = 0.;
      // DIS variables are found in XCM frame, so boost back again
      if (isDIS){
        for (int i = 0; i < 2; ++i) jets[i].transformBy(xcmboost);
      }
      // Note sign change wrt. H1 because photon is in +z direction
      // Jets in HCM so no need to consider orientation.
      for (int i = 0; i < 2; ++i){
        EpPzJets += jets[i].E() - jets[i].pz(); // Sign: + => -
        EmPzJets += jets[i].E() + jets[i].pz(); // Sign: - => +
      }

      // Transform the jets from HCM to LAB frame where eta cut is
      // applied for photoproduction.
      const LorentzTransform hcmboost = kin.boostHCM();
      for (int i = 0; i < 2; ++i) jets[i].transformBy(hcmboost.inverse());
      double etaLabJet1 = dir * jets[0].eta();
      double etaLabJet2 = dir * jets[1].eta();
      double etaMin     = (isDIS) ? -3. : -1.;
      double etaMax     = (isDIS) ? 0. : 2.;
      double eta1       = (isDIS) ? etaXCMJet1 : etaLabJet1;
      double eta2       = (isDIS) ? etaXCMJet2 : etaLabJet2;
      if (!inRange(eta1, etaMin, etaMax)) vetoEvent;
      if (!inRange(eta2, etaMin, etaMax)) vetoEvent;
      ++nVeto5;

      // Pseudorapidity distributions are examined in lab frame:
      double avgEtaJets   = 0.5 * (etaLabJet1 + etaLabJet2);

      // Derive xPom and xGam values from the jet kinematics.
      double zPomJets, xGamJets;
      if (isDIS) {
        zPomJets = (Q2 + M2jets) / (Q2 + M2X);
        xGamJets = EmPzJets / rg.EmPzX(RapidityGap::XCM);
      } else {
        // Boost E_p, E_e to HCM frame
        FourMomentum lep = hcmboost.transform(kin.beamLepton().momentum());
        FourMomentum had = hcmboost.transform(kin.beamHadron().momentum());
        zPomJets = EpPzJets / (2. * xPom * had.E());
        xGamJets = EmPzJets / (2. * y * lep.E());
      }

      // Now fill histograms
      if (isDIS){
        _h_DIS_dsigdzPom     ->fill(zPomJets);
        _h_DIS_dsigdlogXpom  ->fill(log10(xPom));
        _h_DIS_dsigdW        ->fill(W);
        _h_DIS_dsigdQ2       ->fill(Q2);
        _h_DIS_dsigdEtJet1   ->fill(EtJet1);
        _h_DIS_dsigdAvgEta   ->fill(avgEtaJets);
        _h_DIS_dsigdDeltaEta ->fill(deltaEtaJets);
      } else {
        _h_PHO_dsigdzPom     ->fill(zPomJets);
        _h_PHO_dsigdxGam     ->fill(xGamJets);
        _h_PHO_dsigdlogXpom  ->fill(log10(xPom));
        _h_PHO_dsigdW        ->fill(W);
        _h_PHO_dsigdEtJet1   ->fill(EtJet1);
        _h_PHO_dsigdAvgEta   ->fill(avgEtaJets);
        _h_PHO_dsigdDeltaEta ->fill(deltaEtaJets);
        _h_PHO_dsigdMjets    ->fill(sqrt(M2jets));
      }

    }

    // Finalize
    void finalize() {
      // Normalise to cross section
      const double norm = crossSection()/picobarn/sumOfWeights();

      scale( _h_DIS_dsigdzPom    , norm);
      scale( _h_DIS_dsigdlogXpom , norm);
      scale( _h_DIS_dsigdW       , norm);
      scale( _h_DIS_dsigdQ2      , norm);
      scale( _h_DIS_dsigdEtJet1  , norm);
      scale( _h_DIS_dsigdAvgEta  , norm);
      scale( _h_DIS_dsigdDeltaEta, norm);

      scale( _h_PHO_dsigdzPom    , norm);
      scale( _h_PHO_dsigdxGam    , norm);
      scale( _h_PHO_dsigdlogXpom , norm);
      scale( _h_PHO_dsigdW       , norm);
      scale( _h_PHO_dsigdEtJet1  , norm);
      scale( _h_PHO_dsigdAvgEta  , norm);
      scale( _h_PHO_dsigdDeltaEta, norm);
      scale( _h_PHO_dsigdMjets   , norm);

      const double dPHO = nPHO;
      MSG_INFO("H1_2007_I746380");
      MSG_INFO("Cross section = " << crossSection()/picobarn << " pb");
      MSG_INFO("Number of events = " << numEvents() << ", sumW = " << sumOfWeights());
      MSG_INFO("Number of PHO = " << nPHO << ", number of DIS = " << nDIS);
      MSG_INFO("Events passing electron veto   = " << nVeto0 << " (" << nVeto0/dPHO * 100. << "%)" );
      MSG_INFO("Events passing MY              = " << nVeto1 << " (" << nVeto1/dPHO * 100. << "%)" );
      MSG_INFO("Events passing t veto          = " << nVeto2 << " (" << nVeto2/dPHO * 100. << "%)" );
      MSG_INFO("Events passing xPom            = " << nVeto3 << " (" << nVeto3/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet Et veto     = " << nVeto4 << " (" << nVeto4/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet eta veto    = " << nVeto5 << " (" << nVeto5/dPHO * 100. << "%)" );

    }

    //@}


  private:

    /// @name Histograms
    //@{
    // Book histograms from REF data
    Histo1DPtr _h_DIS_dsigdzPom    ;
    Histo1DPtr _h_DIS_dsigdlogXpom ;
    Histo1DPtr _h_DIS_dsigdW       ;
    Histo1DPtr _h_DIS_dsigdQ2      ;
    Histo1DPtr _h_DIS_dsigdEtJet1  ;
    Histo1DPtr _h_DIS_dsigdAvgEta  ;
    Histo1DPtr _h_DIS_dsigdDeltaEta;

    Histo1DPtr _h_PHO_dsigdzPom    ;
    Histo1DPtr _h_PHO_dsigdxGam    ;
    Histo1DPtr _h_PHO_dsigdlogXpom ;
    Histo1DPtr _h_PHO_dsigdW       ;
    Histo1DPtr _h_PHO_dsigdEtJet1  ;
    Histo1DPtr _h_PHO_dsigdAvgEta  ;
    Histo1DPtr _h_PHO_dsigdDeltaEta;
    Histo1DPtr _h_PHO_dsigdMjets   ;
    //@}

    bool isDIS;
    int  nVeto0, nVeto1, nVeto2, nVeto3, nVeto4, nVeto5;
    int nPHO, nDIS;
  };

  DECLARE_RIVET_PLUGIN(H1_2007_I746380);

}

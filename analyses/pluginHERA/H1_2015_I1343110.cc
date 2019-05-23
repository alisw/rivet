// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/DISDiffHadron.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

namespace H1_2015_I1343110_PROJECTIONS {
  
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
      addProjection(DISKinematics(), "DISKIN");
      addProjection(DISFinalState(DISFinalState::HCM), "DISFS");
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

    virtual int compare(const Projection& p) const {
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
      foreach (const Particle& ip, particles) {
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
      foreach (const Particle& jp, pX) {
        momX  += jp.momentum();
        _ePpzX_HCM += jp.E() - jp.pz(); // Sign + => -
        _eMpzX_HCM += jp.E() + jp.pz(); // Sign - => +
      }
      _momX_HCM = momX;
      _pX_HCM   = pX;
      _M2X      = _momX_HCM.mass2();
      
      // Y - side
      FourMomentum momY;
      foreach (const Particle& kp, pY) momY += kp.momentum();
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
      
      foreach (const Particle& jp, pX) {
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

      foreach (const Particle& jp, pY) {
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
      addProjection(RapidityGap(), "RAPGAP");
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
    int compare(const Projection& p) const {
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
  /// Tagged protons & jets found in gamma*p rest frame.
  ///
  /// @author Christine O. Rasmussen
  class H1_2015_I1343110 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(H1_2015_I1343110);

    typedef H1_2015_I1343110_PROJECTIONS::RapidityGap RapidityGap;
    typedef H1_2015_I1343110_PROJECTIONS::BoostedXSystem BoostedXSystem;

    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {

      declare(DISKinematics(), "Kinematics");
      const DISFinalState& disfs = declare(DISFinalState(DISFinalState::HCM), "DISFS");
      const BoostedXSystem& disfsXcm = declare( BoostedXSystem(disfs), "BoostedXFS");
      declare(FastJets(disfsXcm, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 1.0, 
        JetAlg::ALL_MUONS, JetAlg::NO_INVISIBLES, nullptr), "DISFSJets");
      declare(DISDiffHadron(), "Hadron");
      declare(RapidityGap(), "RapidityGap");

      // Book histograms from REF data
      _h_PHO_sig_sqrts     = bookHisto1D(1, 1, 1);
      _h_DIS_sig_sqrts     = bookHisto1D(2, 1, 1);
      _h_PHODIS_sqrts      = bookScatter2D(3, 1, 1);

      _h_DIS_dsigdz        = bookHisto1D(4, 1, 1);
      _h_DIS_dsigdxPom     = bookHisto1D(5, 1, 1);
      _h_DIS_dsigdy        = bookHisto1D(6, 1, 1);
      _h_DIS_dsigdQ2       = bookHisto1D(7, 1, 1);
      _h_DIS_dsigdEtj1     = bookHisto1D(8, 1, 1);
      _h_DIS_dsigdMX       = bookHisto1D(9, 1, 1);
      _h_DIS_dsigdDeltaEta = bookHisto1D(10, 1, 1);
      _h_DIS_dsigdAvgEta   = bookHisto1D(11, 1, 1);

      _h_PHO_dsigdz        = bookHisto1D(12, 1, 1);
      _h_PHO_dsigdxPom     = bookHisto1D(13, 1, 1);
      _h_PHO_dsigdy        = bookHisto1D(14, 1, 1);
      _h_PHO_dsigdxGam     = bookHisto1D(15, 1, 1);
      _h_PHO_dsigdEtj1     = bookHisto1D(16, 1, 1);
      _h_PHO_dsigdMX       = bookHisto1D(17, 1, 1);
      _h_PHO_dsigdDeltaEta = bookHisto1D(18, 1, 1);
      _h_PHO_dsigdAvgEta   = bookHisto1D(19, 1, 1);

      _h_PHODIS_deltaEta   = bookScatter2D(20, 1, 1);
      _h_PHODIS_y          = bookScatter2D(21, 1, 1);
      _h_PHODIS_z          = bookScatter2D(22, 1, 1);
      _h_PHODIS_Etj1       = bookScatter2D(23, 1, 1);
      
      isPHO  = false;
      nVeto1 = 0;
      nVeto2 = 0;
      nVeto3 = 0;
      nVeto4 = 0;
      nVeto5 = 0;
      nVeto6 = 0;
      nPHO   = 0;
      nDIS   = 0;
    }

    // Do the analysis
    void analyze(const Event& event) {

      // Event weight
      const double weight = event.weight(); 
      isPHO  = false;
      
      // Projections - special handling of events where no proton found:
      const RapidityGap&    rg = apply<RapidityGap>(event, "RapidityGap");
      const DISKinematics& kin = apply<DISKinematics>(event, "Kinematics");
      const BoostedXSystem& disfsXcm = apply<BoostedXSystem>( event, "BoostedXFS");
      Particle hadronOut;
      Particle hadronIn;   
      try {
      	const DISDiffHadron& diffhadr = apply<DISDiffHadron>(event, "Hadron");
        hadronOut = diffhadr.out();
        hadronIn  = diffhadr.in();
      } catch (const Error& e){
        vetoEvent;
      }

      // Determine kinematics: H1 has +z = proton direction
      int dir   = kin.orientation();
      double y  = kin.y();
      double Q2 = kin.Q2();

      // Separate into DIS and PHO regimes else veto
      if (Q2 < 2.*GeV2 && inRange(y, 0.2, 0.70)) {
        isPHO = true;
        ++nPHO;
      } else if (inRange(Q2, 4.0*GeV2, 80.*GeV2) && inRange(y, 0.2, 0.7)) {
        isPHO = false;
        ++nDIS;
      } else vetoEvent;
      ++nVeto1;

      // Find diffractive variables as defined in paper. 
      // Note tagged protons in VFPS => smaller allowed xPom range
      // xPom = 1 - E'/E, M2X from hadrons, t = (P-P')^2
      const double M2X  = rg.M2X();
      const double abst = abs(rg.t());
      const double xPom = 1. - hadronOut.energy() / hadronIn.energy();
      
    //cout << "\nhadout=" << hadronOut.energy() << ", hadin=" << hadronIn.energy() << endl;
    //cout << "xPomH1=" << (Q2+M2X) / (y * sqr(sqrtS())) << endl;
    //cout << "|t|=" << abst << ", xPom=" << xPom << endl;
      // Veto if outside allowed region
      if (abst > 0.6 * GeV2)              vetoEvent;
      ++nVeto2;
      if (!inRange(xPom, 0.010, 0.024)) vetoEvent;
      ++nVeto3;

      // Jet selection. Note jets are found in XCM frame, but 
      // eta cut is applied in lab frame! 
      Cut jetcuts = Cuts::Et > 4.* GeV;
      Jets jets   = apply<FastJets>(event, "DISFSJets").jets(jetcuts, cmpMomByEt);
      // Veto if not dijets and if Et_j1 < 5.5
      if (jets.size() < 2)          vetoEvent;
      if (jets[0].Et() < 5.5 * GeV) vetoEvent;
      ++nVeto4;
      // Find Et_jet1 in XCM frame
      double EtJet1 = jets[0].Et() * GeV;
      
    //cout << "gamma*p frame:" << endl;
    //cout << "Et1=" << jets[0].Et() << ", E1=" << jets[0].E() << ", pz1=" << jets[0].pz()  << ", m1=" << jets[0].mass() << endl;
    //cout << "Et2=" << jets[1].Et() << ", E2=" << jets[1].E() << ", pz2=" << jets[1].pz()  << ", m2=" << jets[1].mass() << endl;
      
      // Transform from XCM to HCM
      const LorentzTransform xcmboost = disfsXcm.boost();
      for (int i = 0; i < 2; ++i) jets[i].transformBy(xcmboost.inverse());
      
      // Find mass of jets and EpPz, EmPz of jets in HCM frame.
      FourMomentum momJets = jets[0].momentum() + jets[1].momentum();
      double M2jets        = momJets.mass2();
      double EpPzJets      = 0.;
      double EmPzJets      = 0.;
      // Note sign change wrt. H1 because photon is in +z direction
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
      if (!inRange(etaLabJet1, -1., 2.5)) vetoEvent;
      if (!inRange(etaLabJet2, -1., 2.5)) vetoEvent;
      ++nVeto5;

      // Pseudorapidity distributions are examined in lab frame:      
      double deltaEtaJets = abs(dir * jets[0].eta() - dir * jets[1].eta());
      double avgEtaJets   = 0.5 * (dir * jets[0].eta() + dir * jets[1].eta());
      
      // Evaluate observables
      double zPomJets, xGamJets = 0.;
      if (isPHO){
        zPomJets = EpPzJets / rg.EpPzX(RapidityGap::HCM);
        xGamJets = EmPzJets / rg.EmPzX(RapidityGap::HCM);
        //cout << "xGamJets=" << xGamJets << endl;
      } else { 
        zPomJets = (Q2 + M2jets) / (Q2 + M2X);
      }

    //cout << "lab frame:" << endl;     
    //cout << "Et1=" << jets[0].Et() << ", E1=" << jets[0].E() << ", pz1=" << jets[0].pz() << ", m1=" << jets[0].mass() <<  endl;
    //cout << "Et2=" << jets[1].Et() << ", E2=" << jets[1].E() << ", pz2=" << jets[1].pz() << ", m2=" << jets[1].mass() <<  endl;
    //cout << "EpPzJets=" << EpPzJets << ", EmPzJets=" << EmPzJets << endl;
    //cout << "Et*exp(eta)=" << jets[0].Et()*exp(etaLabJet1) + jets[1].Et()*exp(etaLabJet2) << endl;
    //cout << "Et*exp(-eta)=" << jets[0].Et()*exp(-etaLabJet1) + jets[1].Et()*exp(-etaLabJet2) << endl;
    //cout << "EpPz=" << rg.EpPzX(RapidityGap::HCM) << ", EmPz=" << rg.EmPzX(RapidityGap::HCM) << endl;
    //cout << "2 xPom Ep=" << 2. * xPom * kin.beamHadron().E() << ", 2 y Ee=" << 2. * y * kin.beamLepton().E() << endl;
    //cout << "xGam=" << xGamJets << ", zPom=" << zPomJets << endl;
    //cout << "M12=" << M2jets << ", deltaEta=" << deltaEtaJets << ", avgEta=" << avgEtaJets << endl;
      
      // Veto events with zPom > 0.8
      if (zPomJets > 0.8) vetoEvent;
      ++nVeto6;

      // Now fill histograms
      if (isPHO){
        _h_PHO_sig_sqrts     ->fill(sqrtS()/GeV,   weight);
        _h_PHO_dsigdz        ->fill(zPomJets,      weight);
        _h_PHO_dsigdxPom     ->fill(xPom,          weight);
        _h_PHO_dsigdy        ->fill(y,             weight);
        _h_PHO_dsigdxGam     ->fill(xGamJets,      weight);
        _h_PHO_dsigdEtj1     ->fill(EtJet1,        weight);
        _h_PHO_dsigdMX       ->fill(sqrt(M2X)*GeV, weight);
        _h_PHO_dsigdDeltaEta ->fill(deltaEtaJets,  weight);
        _h_PHO_dsigdAvgEta   ->fill(avgEtaJets,    weight);
      } else {
      	_h_DIS_sig_sqrts     ->fill(sqrtS()/GeV,   weight);
        _h_DIS_dsigdz        ->fill(zPomJets,      weight);
        _h_DIS_dsigdxPom     ->fill(xPom,          weight);
        _h_DIS_dsigdy        ->fill(y,             weight);
        _h_DIS_dsigdQ2       ->fill(Q2,            weight);
        _h_DIS_dsigdEtj1     ->fill(EtJet1,        weight);
        _h_DIS_dsigdMX       ->fill(sqrt(M2X)*GeV, weight);
        _h_DIS_dsigdDeltaEta ->fill(deltaEtaJets,  weight);
        _h_DIS_dsigdAvgEta   ->fill(avgEtaJets,    weight);
      }                    
      
    }

    // Finalize
    void finalize() {
      // Normalise to cross section
      // Remember to manually scale the cross section afterwards with
      // the number of rejected events.
      const double norm = crossSection()/picobarn/sumOfWeights();
      
      scale(_h_PHO_sig_sqrts,     norm); 
      scale(_h_PHO_dsigdz,        norm); 
      scale(_h_PHO_dsigdxPom,     norm); 
      scale(_h_PHO_dsigdy,        norm); 
      scale(_h_PHO_dsigdxGam,     norm); 
      scale(_h_PHO_dsigdEtj1,     norm); 
      scale(_h_PHO_dsigdMX,       norm); 
      scale(_h_PHO_dsigdDeltaEta, norm); 
      scale(_h_PHO_dsigdAvgEta,   norm); 

      scale(_h_DIS_sig_sqrts,     norm); 
      scale(_h_DIS_dsigdz,        norm); 
      scale(_h_DIS_dsigdxPom,     norm); 
      scale(_h_DIS_dsigdy,        norm); 
      scale(_h_DIS_dsigdQ2,       norm); 
      scale(_h_DIS_dsigdEtj1,     norm); 
      scale(_h_DIS_dsigdMX,       norm); 
      scale(_h_DIS_dsigdDeltaEta, norm); 
      scale(_h_DIS_dsigdAvgEta,   norm); 
      
      if (_h_DIS_sig_sqrts->numEntries() != 0)
        divide(_h_PHO_sig_sqrts, _h_DIS_sig_sqrts, _h_PHODIS_sqrts); 
      if (_h_DIS_dsigdDeltaEta->numEntries() != 0)
        divide(_h_PHO_dsigdDeltaEta, _h_DIS_dsigdDeltaEta, _h_PHODIS_deltaEta); 
      if (_h_DIS_dsigdy->numEntries() != 0)
        divide(_h_PHO_dsigdy, _h_DIS_dsigdy, _h_PHODIS_y); 
      if (_h_DIS_dsigdz->numEntries() != 0)
        divide(_h_PHO_dsigdz, _h_DIS_dsigdz, _h_PHODIS_z); 
      if (_h_DIS_dsigdEtj1->numEntries() != 0)
        divide(_h_PHO_dsigdEtj1, _h_DIS_dsigdEtj1, _h_PHODIS_Etj1); 

      const double dPHO = nPHO;
      MSG_INFO("H1_2015_I1343110");
      MSG_INFO("Cross section = " << crossSection()/picobarn << " pb");
      MSG_INFO("Number of events = " << numEvents() << ", sumW = " << sumOfWeights());
      MSG_INFO("Number of PHO = " << nPHO << ", number of DIS = " << nDIS);
      MSG_INFO("Events passing electron veto   = " << nVeto1 << " (" << nVeto1/dPHO * 100. << "%)" );
      MSG_INFO("Events passing t veto          = " << nVeto2 << " (" << nVeto2/dPHO * 100. << "%)" );
      MSG_INFO("Events passing xPom            = " << nVeto3 << " (" << nVeto3/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet Et   veto   = " << nVeto4 << " (" << nVeto4/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet eta veto    = " << nVeto5 << " (" << nVeto5/dPHO * 100. << "%)" );
      MSG_INFO("Events passing zPom veto       = " << nVeto6 << " (" << nVeto6/dPHO * 100. << "%)" );

    }

    //@}


  private:

    /// @name Histograms
    //@{
    // Book histograms from REF data
    Histo1DPtr _h_PHO_sig_sqrts;
    Histo1DPtr _h_DIS_sig_sqrts;
    Scatter2DPtr _h_PHODIS_sqrts;
    
    Histo1DPtr _h_DIS_dsigdz;
    Histo1DPtr _h_DIS_dsigdxPom;
    Histo1DPtr _h_DIS_dsigdy;
    Histo1DPtr _h_DIS_dsigdQ2;
    Histo1DPtr _h_DIS_dsigdEtj1;
    Histo1DPtr _h_DIS_dsigdMX;
    Histo1DPtr _h_DIS_dsigdDeltaEta;
    Histo1DPtr _h_DIS_dsigdAvgEta;

    Histo1DPtr _h_PHO_dsigdz;
    Histo1DPtr _h_PHO_dsigdxPom;
    Histo1DPtr _h_PHO_dsigdy;
    Histo1DPtr _h_PHO_dsigdxGam;
    Histo1DPtr _h_PHO_dsigdEtj1;
    Histo1DPtr _h_PHO_dsigdMX;
    Histo1DPtr _h_PHO_dsigdDeltaEta;
    Histo1DPtr _h_PHO_dsigdAvgEta;

    Scatter2DPtr _h_PHODIS_deltaEta;
    Scatter2DPtr _h_PHODIS_y;
    Scatter2DPtr _h_PHODIS_z;
    Scatter2DPtr _h_PHODIS_Etj1;
    //@}

    bool isPHO;
    int  nVeto1, nVeto2, nVeto3, nVeto4, nVeto5, nVeto6;
    int  nPHO, nDIS;
  };

  DECLARE_RIVET_PLUGIN(H1_2015_I1343110);

}

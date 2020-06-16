// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  /// @brief D0 topological distributions of 3- and 4-jet events.
  class D0_1996_S3214044 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_1996_S3214044() : Analysis("D0_1996_S3214044")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      const FinalState fs;
      declare(fs, "FS");
      /// @todo Use correct jet algorithm --- tried FJ3 D0RunICone but does
      // not look as good as the Run2 cone alg used here
      declare(FastJets(fs, FastJets::D0ILCONE, 0.7), "ConeJets");

      book(_h_3j_x3 ,1, 1, 1);
      book(_h_3j_x5 ,2, 1, 1);
      book(_h_3j_costheta3 ,3, 1, 1);
      book(_h_3j_psi ,4, 1, 1);
      book(_h_3j_mu34 ,5, 1, 1);
      book(_h_3j_mu35 ,6, 1, 1);
      book(_h_3j_mu45 ,7, 1, 1);

      book(_h_4j_x3 ,8, 1, 1);
      book(_h_4j_x4 ,9, 1, 1);
      book(_h_4j_x5 ,10, 1, 1);
      book(_h_4j_x6 ,11, 1, 1);
      book(_h_4j_costheta3 ,12, 1, 1);
      book(_h_4j_costheta4 ,13, 1, 1);
      book(_h_4j_costheta5 ,14, 1, 1);
      book(_h_4j_costheta6 ,15, 1, 1);
      book(_h_4j_cosomega34 ,16, 1, 1);
      book(_h_4j_cosomega35 ,17, 1, 1);
      book(_h_4j_cosomega36 ,18, 1, 1);
      book(_h_4j_cosomega45 ,19, 1, 1);
      book(_h_4j_cosomega46 ,20, 1, 1);
      book(_h_4j_cosomega56 ,21, 1, 1);
      book(_h_4j_mu34 ,22, 1, 1);
      book(_h_4j_mu35 ,23, 1, 1);
      book(_h_4j_mu36 ,24, 1, 1);
      book(_h_4j_mu45 ,25, 1, 1);
      book(_h_4j_mu46 ,26, 1, 1);
      book(_h_4j_mu56 ,27, 1, 1);
      book(_h_4j_theta_BZ ,28, 1, 1);
      book(_h_4j_costheta_NR ,29, 1, 1);

    }


    void analyze(const Event& event) {
      Jets jets_in = apply<FastJets>(event, "ConeJets")
        .jets(Cuts::Et > 20*GeV && Cuts::abseta < 3, cmpMomByEt);

      Jets jets_isolated;
      for (size_t i = 0; i < jets_in.size(); ++i) {
        bool isolated = true;
        for (size_t j = 0; j < jets_in.size(); ++j) {
          if (i != j && deltaR(jets_in[i], jets_in[j]) < 1.4) {
            isolated = false;
            break;
          }
        }
        if (isolated) jets_isolated.push_back(jets_in[i]);
      }

      if (jets_isolated.size() == 0 || jets_isolated[0].Et() < 60.0*GeV) vetoEvent;

      if (jets_isolated.size() > 2) _threeJetAnalysis(jets_isolated);
      if (jets_isolated.size() > 3) _fourJetAnalysis(jets_isolated);
    }


    void finalize() {
      normalize(_h_3j_x3, 1.0);
      normalize(_h_3j_x5, 1.0);
      normalize(_h_3j_costheta3, 1.0);
      normalize(_h_3j_psi, 1.0);
      normalize(_h_3j_mu34, 1.0);
      normalize(_h_3j_mu35, 1.0);
      normalize(_h_3j_mu45, 1.0);
      normalize(_h_4j_x3, 1.0);
      normalize(_h_4j_x4, 1.0);
      normalize(_h_4j_x5, 1.0);
      normalize(_h_4j_x6, 1.0);
      normalize(_h_4j_costheta3, 1.0);
      normalize(_h_4j_costheta4, 1.0);
      normalize(_h_4j_costheta5, 1.0);
      normalize(_h_4j_costheta6, 1.0);
      normalize(_h_4j_cosomega34, 1.0);
      normalize(_h_4j_cosomega35, 1.0);
      normalize(_h_4j_cosomega36, 1.0);
      normalize(_h_4j_cosomega45, 1.0);
      normalize(_h_4j_cosomega46, 1.0);
      normalize(_h_4j_cosomega56, 1.0);
      normalize(_h_4j_mu34, 1.0);
      normalize(_h_4j_mu35, 1.0);
      normalize(_h_4j_mu36, 1.0);
      normalize(_h_4j_mu45, 1.0);
      normalize(_h_4j_mu46, 1.0);
      normalize(_h_4j_mu56, 1.0);
      normalize(_h_4j_theta_BZ, 1.0);
      normalize(_h_4j_costheta_NR, 1.0);
    }

    //@}


  private:

    /// @name Helper functions
    //@{

    void _threeJetAnalysis(const Jets& jets) {
      // >=3 jet events
      FourMomentum jjj(jets[0].momentum()+jets[1].momentum()+jets[2].momentum());
      const double sqrts = _safeMass(jjj);
      if (sqrts<200*GeV) {
        return;
      }

      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(jjj.betaVec());
      vector<FourMomentum> jets_boosted;
      for (Jet jet : jets) {
        jets_boosted.push_back(cms_boost.transform(jet.momentum()));
      }
      std::sort(jets_boosted.begin(), jets_boosted.end(), FourMomentum::byEDescending());
      FourMomentum p3(jets_boosted[0]);
      FourMomentum p4(jets_boosted[1]);
      FourMomentum p5(jets_boosted[2]);

      Vector3 beam1(0.0, 0.0, 1.0);
      Vector3 p1xp3 = beam1.cross(p3.p3());
      Vector3 p4xp5 = p4.p3().cross(p5.p3());
      const double cospsi = p1xp3.dot(p4xp5)/p1xp3.mod()/p4xp5.mod();

      _h_3j_x3->fill(2.0*p3.E()/sqrts);
      _h_3j_x5->fill(2.0*p5.E()/sqrts);
      _h_3j_costheta3->fill(fabs(cos(p3.theta())));
      _h_3j_psi->fill(acos(cospsi)/degree);
      _h_3j_mu34->fill(_safeMass(FourMomentum(p3+p4))/sqrts);
      _h_3j_mu35->fill(_safeMass(FourMomentum(p3+p5))/sqrts);
      _h_3j_mu45->fill(_safeMass(FourMomentum(p4+p5))/sqrts);
    }


    void _fourJetAnalysis(const Jets& jets) {
      // >=4 jet events
      FourMomentum jjjj(jets[0].momentum() + jets[1].momentum() + jets[2].momentum()+ jets[3].momentum());
      const double sqrts = _safeMass(jjjj);
      if (sqrts < 200*GeV) return;

      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(jjjj.betaVec());
      vector<FourMomentum> jets_boosted;
      for (Jet jet : jets) {
        jets_boosted.push_back(cms_boost.transform(jet.momentum()));
      }
      sort(jets_boosted.begin(), jets_boosted.end(), FourMomentum::byEDescending());
      FourMomentum p3(jets_boosted[0]);
      FourMomentum p4(jets_boosted[1]);
      FourMomentum p5(jets_boosted[2]);
      FourMomentum p6(jets_boosted[3]);

      Vector3 p3xp4 = p3.p3().cross(p4.p3());
      Vector3 p5xp6 = p5.p3().cross(p6.p3());
      const double costheta_BZ = p3xp4.dot(p5xp6)/p3xp4.mod()/p5xp6.mod();
      const double costheta_NR = (p3.p3()-p4.p3()).dot(p5.p3()-p6.p3())/
        (p3.p3()-p4.p3()).mod()/(p5.p3()-p6.p3()).mod();

      _h_4j_x3->fill(2.0*p3.E()/sqrts);
      _h_4j_x4->fill(2.0*p4.E()/sqrts);
      _h_4j_x5->fill(2.0*p5.E()/sqrts);
      _h_4j_x6->fill(2.0*p6.E()/sqrts);
      _h_4j_costheta3->fill(fabs(cos(p3.theta())));
      _h_4j_costheta4->fill(fabs(cos(p4.theta())));
      _h_4j_costheta5->fill(fabs(cos(p5.theta())));
      _h_4j_costheta6->fill(fabs(cos(p6.theta())));
      _h_4j_cosomega34->fill(cos(p3.angle(p4)));
      _h_4j_cosomega35->fill(cos(p3.angle(p5)));
      _h_4j_cosomega36->fill(cos(p3.angle(p6)));
      _h_4j_cosomega45->fill(cos(p4.angle(p5)));
      _h_4j_cosomega46->fill(cos(p4.angle(p6)));
      _h_4j_cosomega56->fill(cos(p5.angle(p6)));
      _h_4j_mu34->fill(_safeMass(FourMomentum(p3+p4))/sqrts);
      _h_4j_mu35->fill(_safeMass(FourMomentum(p3+p5))/sqrts);
      _h_4j_mu36->fill(_safeMass(FourMomentum(p3+p6))/sqrts);
      _h_4j_mu45->fill(_safeMass(FourMomentum(p4+p5))/sqrts);
      _h_4j_mu46->fill(_safeMass(FourMomentum(p4+p6))/sqrts);
      _h_4j_mu56->fill(_safeMass(FourMomentum(p5+p6))/sqrts);
      _h_4j_theta_BZ->fill(acos(fabs(costheta_BZ))/degree);
      _h_4j_costheta_NR->fill(fabs(costheta_NR));

    }

    double _safeMass(const FourMomentum& p) {
      double mass2=p.mass2();
      if (mass2>0.0) return sqrt(mass2);
      else if (mass2<-1.0e-5) {
        MSG_WARNING("m2 = " << m2 << ". Assuming m2=0.");
        return 0.0;
      }
      else return 0.0;
    }

  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_3j_x3;
    Histo1DPtr _h_3j_x5;
    Histo1DPtr _h_3j_costheta3;
    Histo1DPtr _h_3j_psi;
    Histo1DPtr _h_3j_mu34;
    Histo1DPtr _h_3j_mu35;
    Histo1DPtr _h_3j_mu45;

    Histo1DPtr _h_4j_x3;
    Histo1DPtr _h_4j_x4;
    Histo1DPtr _h_4j_x5;
    Histo1DPtr _h_4j_x6;
    Histo1DPtr _h_4j_costheta3;
    Histo1DPtr _h_4j_costheta4;
    Histo1DPtr _h_4j_costheta5;
    Histo1DPtr _h_4j_costheta6;
    Histo1DPtr _h_4j_cosomega34;
    Histo1DPtr _h_4j_cosomega35;
    Histo1DPtr _h_4j_cosomega36;
    Histo1DPtr _h_4j_cosomega45;
    Histo1DPtr _h_4j_cosomega46;
    Histo1DPtr _h_4j_cosomega56;
    Histo1DPtr _h_4j_mu34;
    Histo1DPtr _h_4j_mu35;
    Histo1DPtr _h_4j_mu36;
    Histo1DPtr _h_4j_mu45;
    Histo1DPtr _h_4j_mu46;
    Histo1DPtr _h_4j_mu56;
    Histo1DPtr _h_4j_theta_BZ;
    Histo1DPtr _h_4j_costheta_NR;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_1996_S3214044);

}

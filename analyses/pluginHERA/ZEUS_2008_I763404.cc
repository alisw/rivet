// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISDiffHadron.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ZEUS dijet photoproduction study used in the ZEUS jets PDF fit
  ///
  /// This class is a reproduction of the HZTool routine for the ZEUS
  /// dijet photoproduction paper which was used in the ZEUS jets PDF fit.
  ///
  /// @author Ilkka Helenius
  class ZEUS_2008_I763404 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ZEUS_2008_I763404);

    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {

      /// @todo Acceptance
      FinalState fs;
      // Final state particles with central tracking detector.
      declare(FastJets(fs, FastJets::KT, 1.0), "Jets");

      // Projections
      declare(DISKinematics(), "Kinematics");
      declare(DISDiffHadron(), "Hadron");

      _h_dsigma_all[0] = bookHisto1D(1, 1, 1);
      _h_dsigma_all[1] = bookHisto1D(2, 1, 1);
      _h_dsigma_all[2] = bookHisto1D(3, 1, 1);
      _h_dsigma_all[3] = bookHisto1D(4, 1, 1);
      _h_dsigma_all[4] = bookHisto1D(5, 1, 1);
      _h_dsigma_all[5] = bookHisto1D(6, 1, 1);
      _h_xgamma = bookHisto1D(7, 1, 1);
      _h_dsigma_xgamma[0][0] = bookHisto1D(8, 1, 1);
      _h_dsigma_xgamma[0][1] = bookHisto1D(9, 1, 1);
      _h_dsigma_xgamma[0][2] = bookHisto1D(10, 1, 1);
      _h_dsigma_xgamma[0][3] = bookHisto1D(11, 1, 1);
      _h_dsigma_xgamma[1][0] = bookHisto1D(12, 1, 1);
      _h_dsigma_xgamma[1][1] = bookHisto1D(13, 1, 1);
      _h_dsigma_xgamma[1][2] = bookHisto1D(14, 1, 1);
      _h_dsigma_xgamma[1][3] = bookHisto1D(15, 1, 1);

      nVeto0 = 0;
      nVeto1 = 0;
      nVeto2 = 0;
      nVeto3 = 0;
      nVeto4 = 0;
    }

    // Do the analysis
    void analyze(const Event& event) {

      // Derive the DIS kinematics.
      const DISKinematics& kin   = apply<DISKinematics>(event, "Kinematics");

      // Derive the diffractive kinematics (should be used for diffractive only).
      Particle hadronOut;
      Particle hadronIn;
      try {
      	const DISDiffHadron & diffhadr = apply<DISDiffHadron>(event, "Hadron");
        hadronOut = diffhadr.out();
        hadronIn  = diffhadr.in();
      } catch (const Error& e){
        vetoEvent;
      }

      // Determine event orientation, since coord system is for +z = proton direction
      const int orientation = kin.orientation();

      // Calculate the photon 4-momentum from the incoming and outgoing lepton.
      const FourMomentum qleptonIn  = kin.beamLepton().momentum();
      const FourMomentum qleptonOut = kin.scatteredLepton().momentum();
      const FourMomentum qphoton    = qleptonIn - qleptonOut;

      // Calculate the Pomeron 4-momentum from the incoming and outgoing hadron
      const FourMomentum pHadOut  = hadronOut.momentum();
      const FourMomentum pHadIn   = hadronIn.momentum();
      const FourMomentum pPomeron = pHadIn - pHadOut;

      // Q2 and inelasticity cuts
      if (kin.Q2() > 1*GeV2) vetoEvent;
      ++nVeto0;
      if (!inRange(kin.y(), 0.2, 0.85)) vetoEvent;
      ++nVeto1;

      // Jet selection and veto.
      const Jets jets = apply<FastJets>(event, "Jets") \
        .jets(Cuts::Et > 6.5*GeV && Cuts::etaIn(-1.5*orientation, 1.5*orientation), cmpMomByEt);
      MSG_DEBUG("Jet multiplicity = " << jets.size());
      if (jets.size() < 2) vetoEvent;
      ++nVeto2;
      const Jet& j1 = jets[0];
      const Jet& j2 = jets[1];
      if (j1.Et() < 7.5*GeV) vetoEvent;
      ++nVeto3;

      // Veto on x_Pomeron.
      const double xPom = ( pPomeron * qphoton ) / (pHadIn * qphoton);
      if (xPom > 0.025) vetoEvent;
      ++nVeto4;

      // Computation of event-level variables.
      const double eta1 = orientation*j1.eta(), eta2 = orientation*j2.eta();
      const double xyobs = (j1.Et() * exp(-eta1) + j2.Et() * exp(-eta2)) / (2*kin.y()*kin.beamLepton().E());
      const size_t i_xyobs = (xyobs < 0.75) ? 1 : 0;
      const double zPobs = (j1.Et() * exp(eta1) + j2.Et() * exp(eta2)) / (2*xPom*kin.beamHadron().E());
      const double M_X = sqrt( (pPomeron + qphoton).mass2() );

      // Fill histograms
      const double weight = event.weight();

      _h_dsigma_all[0]->fill(kin.y(), weight);
      _h_dsigma_all[1]->fill(M_X, weight);
      _h_dsigma_all[2]->fill(xPom, weight);
      _h_dsigma_all[3]->fill(zPobs, weight);
      _h_dsigma_all[4]->fill(j1.Et(), weight);
      _h_dsigma_all[5]->fill(eta1, weight);

      _h_xgamma->fill(xyobs, weight);

      _h_dsigma_xgamma[i_xyobs][0]->fill(kin.y(), weight);
      _h_dsigma_xgamma[i_xyobs][1]->fill(M_X, weight);
      _h_dsigma_xgamma[i_xyobs][2]->fill(xPom, weight);
      _h_dsigma_xgamma[i_xyobs][3]->fill(zPobs, weight);

    }

    // Finalize
    void finalize() {
      const double norm = crossSection()/picobarn/sumOfWeights();
      scale(_h_xgamma, norm);
      for (auto& h : _h_dsigma_all) scale(h, norm);
      for (auto& h : _h_dsigma_xgamma[0]) scale(h, norm);
      for (auto& h : _h_dsigma_xgamma[1]) scale(h, norm);

      // Cross section in nb for these observables.
      scale(_h_dsigma_all[2], 1e-3);
      scale(_h_dsigma_xgamma[0][2], 1e-3);
      scale(_h_dsigma_xgamma[1][2], 1e-3);

      double dPHO = nVeto1;

      MSG_INFO("ZEUS_2008_I763403");
      MSG_INFO("Cross section = " << crossSection()/picobarn);
      MSG_INFO("Number of events = " << numEvents() << ", sumW = " << sumOfWeights());
      MSG_INFO("Events passing electron veto1= " << nVeto0 << " (" << nVeto0/dPHO * 100. << "%)" );
      MSG_INFO("Events passing electron veto2= " << nVeto1 << " (" << nVeto1/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet size veto = " << nVeto2 << " (" << nVeto2/dPHO * 100. << "%)" );
      MSG_INFO("Events passing jet Et veto   = " << nVeto3 << " (" << nVeto3/dPHO * 100. << "%)" );
      MSG_INFO("Events passing xPom veto     = " << nVeto4 << " (" << nVeto4/dPHO * 100. << "%)" );

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_dsigma_all[6];
    Histo1DPtr _h_xgamma;
    Histo1DPtr _h_dsigma_xgamma[2][4];
    //@}

    int nVeto0, nVeto1, nVeto2, nVeto3, nVeto4;
  };

  DECLARE_RIVET_PLUGIN(ZEUS_2008_I763404);

}

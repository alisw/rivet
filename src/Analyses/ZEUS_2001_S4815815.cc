// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ZEUS dijet photoproduction study used in the ZEUS jets PDF fit
  ///
  /// This class is a reproduction of the HZTool routine for the ZEUS
  /// dijet photoproduction paper which was used in the ZEUS jets PDF fit.
  ///
  /// @note Cleaning cuts on event pT/sqrt(Et) and y_e are not needed in MC analysis.
  ///
  /// @author Andy Buckley
  class ZEUS_2001_S4815815 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ZEUS_2001_S4815815);

    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {

      /// @todo Acceptance
      FinalState fs;
      declare(FastJets(fs, FastJets::KT, 1.0), "Jets"); //< R=1 checked with Matt Wing

      /// @todo Dress the lepton?
      IdentifiedFinalState positrons(fs, PID::POSITRON);
      declare(positrons, "Positrons");

      // Table 1
      _h_costh[0] = bookHisto1D(1, 1, 1);
      _h_costh[1] = bookHisto1D(1, 1, 2);
      // Table 2
      _h_etjet1[1][0] = bookHisto1D(2, 1, 1);
      _h_etjet1[1][1] = bookHisto1D(3, 1, 1);
      _h_etjet1[1][2] = bookHisto1D(4, 1, 1);
      _h_etjet1[1][3] = bookHisto1D(5, 1, 1);
      _h_etjet1[1][4] = bookHisto1D(6, 1, 1);
      _h_etjet1[1][5] = bookHisto1D(7, 1, 1);
      // Table 3
      _h_etjet1[0][0] = bookHisto1D(8, 1, 1);
      _h_etjet1[0][1] = bookHisto1D(9, 1, 1);
      _h_etjet1[0][2] = bookHisto1D(10, 1, 1);
      _h_etjet1[0][3] = bookHisto1D(11, 1, 1);
      _h_etjet1[0][4] = bookHisto1D(12, 1, 1);
      _h_etjet1[0][5] = bookHisto1D(13, 1, 1);
      // Table 4
      _h_etajet2[1][0] = bookHisto1D(14, 1, 1);
      _h_etajet2[1][1] = bookHisto1D(15, 1, 1);
      _h_etajet2[1][2] = bookHisto1D(16, 1, 1);
      // Table 5
      _h_etajet2[0][0] = bookHisto1D(17, 1, 1);
      _h_etajet2[0][1] = bookHisto1D(18, 1, 1);
      _h_etajet2[0][2] = bookHisto1D(19, 1, 1);
      // Table 6
      _h_xobsy[0] = bookHisto1D(20, 1, 1);
      _h_xobsy[1] = bookHisto1D(21, 1, 1);
      _h_xobsy[2] = bookHisto1D(22, 1, 1);
      _h_xobsy[3] = bookHisto1D(23, 1, 1);
    }


    // Do the analysis
    void analyze(const Event& event) {

      // Determine event orientation, since coord system is for +z = proton direction
      const ParticlePair bs = event.beams();
      if (bs.first.pid() != PID::POSITRON && bs.second.pid() != PID::POSITRON) vetoEvent;
      const Particle& bpositron = (bs.first.pid() == PID::POSITRON ? bs.first : bs.second);
      if (bs.first.pid() != PID::PROTON && bs.second.pid() != PID::PROTON) vetoEvent;
      const Particle& bproton = (bs.first.pid() == PID::PROTON) ? bs.first : bs.second;
      const int orientation = sign(bproton.momentum().pz());
      MSG_DEBUG("Beam proton = " << bproton.mom() << " GeV => orientation = " << orientation);

      // Jet selection
      const Jets jets = apply<FastJets>(event, "Jets") \
        .jets(Cuts::pT > 11*GeV && Cuts::etaIn(-1*orientation, 2.4*orientation), cmpMomByEt);
      MSG_DEBUG("Jet multiplicity = " << jets.size());
      if (jets.size() < 2) vetoEvent;
      const Jet& j1 = jets[0];
      const Jet& j2 = jets[1];
      if (j1.pT() < 14*GeV) vetoEvent;

      // eta and cos(theta*) computation
      const double eta1 = orientation*j1.eta(), eta2 = orientation*j2.eta();
      const double etabar = (eta1 + eta2)/2;
      const double etadiff = eta1 - eta2;
      const double costhetastar = tanh(etadiff/2);

      // Get the scattered positron
      const Particles positrons = apply<FinalState>(event, "Positrons").particlesByPt();
      if (positrons.empty()) vetoEvent;
      const Particle& positron = positrons.front();

      // Calculate the photon 4-vector
      const FourMomentum qphoton = positron.mom() - bpositron.mom();

      // Computation and cut on inelasticity
      const double inelasticity = dot(bproton.mom(), qphoton) / dot(bproton.mom(), bpositron.mom());
      if (!inRange(inelasticity, 0.2, 0.85)) vetoEvent;

      // Computation of x_y^obs
      // (I assume Ee is the lab frame positron momentum, not in proton rest frame cf. the ambiguous phrase in the paper)
      const double xyobs = (j1.Et() * exp(-eta1) + j2.Et() * exp(-eta2)) / (2*inelasticity*bpositron.E());
      const size_t i_xyobs = (xyobs < 0.75) ? 0 : 1;

      // Fill histograms
      const double weight = event.weight();
      // T1
      if ((j1.mom()+j2.mom()).mass() > 42*GeV && inRange(etabar, 0.1, 0.3))
        _h_costh[i_xyobs]->fill(abs(costhetastar), weight);
      // T2, T3
      if (inRange(eta1, -1, 0) && inRange(eta2, -1, 0))
        _h_etjet1[i_xyobs][0]->fill(j1.Et()/GeV, weight);
      else if (inRange(eta1, 0, 1) && inRange(eta2, -1, 0))
        _h_etjet1[i_xyobs][1]->fill(j1.Et()/GeV, weight);
      else if (inRange(eta1, 0, 1) && inRange(eta2, 0, 1))
        _h_etjet1[i_xyobs][2]->fill(j1.Et()/GeV, weight);
      else if (inRange(eta1, 1, 2.4) && inRange(eta2, -1, 0))
        _h_etjet1[i_xyobs][3]->fill(j1.Et()/GeV, weight);
      else if (inRange(eta1, 1, 2.4) && inRange(eta2, 0, 1))
        _h_etjet1[i_xyobs][4]->fill(j1.Et()/GeV, weight);
      else if (inRange(eta1, 1, 2.4) && inRange(eta2, 1, 2.4))
        _h_etjet1[i_xyobs][5]->fill(j1.Et()/GeV, weight);
      // T4, T5
      if (inRange(eta1, -1, 0))
        _h_etajet2[i_xyobs][0]->fill(eta2, weight);
      else if (inRange(eta1, 0, 1))
        _h_etajet2[i_xyobs][1]->fill(eta2, weight);
      else if (inRange(eta1, 1, 2.4))
        _h_etajet2[i_xyobs][2]->fill(eta2, weight);
      // T6
      if (inRange(j1.Et()/GeV, 14, 17))
        _h_xobsy[0]->fill(xyobs, weight);
      else if (inRange(j1.Et()/GeV, 17, 25))
        _h_xobsy[1]->fill(xyobs, weight);
      else if (inRange(j1.Et()/GeV, 25, 35))
        _h_xobsy[2]->fill(xyobs, weight);
      else if (inRange(j1.Et()/GeV, 35, 90))
        _h_xobsy[3]->fill(xyobs, weight);
    }


    // Finalize
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();
      for (size_t ix = 0; ix < 2; ++ix) {
        scale(_h_costh[ix], sf);
        for (auto& h : _h_etjet1[ix]) scale(h, sf);
        for (auto& h : _h_etajet2[ix]) scale(h, sf);
      }
      for (auto& h : _h_xobsy) scale(h, sf);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_costh[2], _h_etjet1[2][6], _h_etajet2[2][3], _h_xobsy[4];
    //@}

  };


  DECLARE_RIVET_PLUGIN(ZEUS_2001_S4815815);

}

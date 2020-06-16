// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  namespace {

    /// @name Jet angle calculator functions
    //@{
    /// @todo Move to utils? (taken from DELPHI_2003)

    /// @todo Use Jet or FourMomentum interface rather than PseudoJet
    /// @todo Move to utils?
    double calc_BZ(const vector<fastjet::PseudoJet>& jets) {
      assert(jets.size() == 4);
      Vector3 p12 = cross( momentum3(jets[0]), momentum3(jets[1]));
      Vector3 p34 = cross( momentum3(jets[2]), momentum3(jets[3]));
      return dot(p12,p34) / (p12.mod()*p34.mod());
    }


    /// @todo Use Jet or FourMomentum interface rather than PseudoJet
    /// @todo Move to utils?
    double calc_KSW(const vector<fastjet::PseudoJet>& jets) {
      assert(jets.size() == 4);
      Vector3 p13 = cross( momentum3(jets[0]), momentum3(jets[2]));
      Vector3 p24 = cross( momentum3(jets[1]), momentum3(jets[3]));
      Vector3 p14 = cross( momentum3(jets[0]), momentum3(jets[3]));
      Vector3 p23 = cross( momentum3(jets[1]), momentum3(jets[2]));
      return cos (0.5*( acos (dot(p14,p23) / (p14.mod()*p23.mod())) +
                        acos (dot(p13,p24) / (p13.mod()*p24.mod())) ));
    }


    /// @todo Use Jet or FourMomentum interface rather than PseudoJet
    /// @todo Move to utils?
    double calc_NR(const vector<fastjet::PseudoJet>& jets) {
      assert(jets.size() == 4);
      Vector3 p12 = momentum3(jets[0]) - momentum3(jets[1]);
      Vector3 p34 = momentum3(jets[2]) - momentum3(jets[3]);
      return dot(p12,p34) / (p12.mod()*p34.mod());
    }

    /// @todo Use Jet or FourMomentum interface rather than PseudoJet
    /// @todo Move to utils?
    double calc_ALPHA34(const vector<fastjet::PseudoJet>& jets) {
      assert(jets.size() == 4);
      Vector3 p3 = momentum3(jets[2]);
      Vector3 p4 = momentum3(jets[3]);
      return dot(p3,p4) / (p3.mod()*p4.mod());
    }

    //@}

  }


  class OPAL_2001_S4553896 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    OPAL_2001_S4553896()
      : Analysis("OPAL_2001_S4553896")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      const FinalState fs;
      declare(fs, "FS");
      declare(FastJets(fs, FastJets::DURHAM, 0.7), "Jets");


      /// @todo Book histograms here, e.g.:
      book(_h_BZ      ,3, 1, 1);
      book(_h_KSW     ,4, 1, 1);
      book(_h_NR      ,5, 1, 1);
      book(_h_ALPHA34 ,6, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
 
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (apply<FinalState>(event, "FS").particles().size() < 2) {
        vetoEvent;
      }

      const FastJets& fastjets = apply<FastJets>(event, "Jets");
      if (fastjets.clusterSeq()) {
        vector<fastjet::PseudoJet> jets;
        for (const fastjet::PseudoJet& jet :
                 fastjet::sorted_by_E(fastjets.clusterSeq()->exclusive_jets_ycut(0.008))) {
          if (jet.E()>3.0*GeV) jets.push_back(jet);
        }
        if (jets.size() == 4) {
          // Prevent nan-fill due to division by zero in calc_BZ
          double bz = fabs(calc_BZ(jets));
          if (!std::isnan(bz)) _h_BZ->fill(bz);
          _h_KSW->fill(calc_KSW(jets));
          _h_NR->fill(fabs(calc_NR(jets)));
          _h_ALPHA34->fill(calc_ALPHA34(jets));
        }
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// Normalise, scale and otherwise manipulate histograms here
      normalize(_h_BZ);
      normalize(_h_KSW);
      normalize(_h_NR);
      normalize(_h_ALPHA34);

    }

    //@}


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_BZ;
    Histo1DPtr _h_KSW;
    Histo1DPtr _h_NR;
    Histo1DPtr _h_ALPHA34;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2001_S4553896);

}

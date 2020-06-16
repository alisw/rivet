// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Jet rates in \f$ e^+ e^- \f$ at OPAL and JADE
  /// @author Frank Siegert
  class JADE_OPAL_2000_S4300807 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_OPAL_2000_S4300807)


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      const FinalState fs;
      declare(fs, "FS");
      FastJets jadeJets = FastJets(fs, FastJets::JADE, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      FastJets durhamJets = FastJets(fs, FastJets::DURHAM, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jadeJets, "JadeJets");
      declare(durhamJets, "DurhamJets");

      // Histos
      int offset = 0;
      switch (int(sqrtS()/GeV + 0.5)) {
      case 35: offset = 7; break;
      case 44: offset = 8; break;
      case 91: offset = 9; break;
      case 133: offset = 10; break;
      case 161: offset = 11; break;
      case 172: offset = 12; break;
      case 183: offset = 13; break;
      case 189: offset = 14; break;
      default: break;
      }
      for (size_t i = 0; i < 5; ++i) {
        book(_h_R_Jade[i] ,offset, 1, i+1);
        book(_h_R_Durham[i] ,offset+9, 1, i+1);
        if (i < 4) book(_h_y_Durham[i] ,offset+17, 1, i+1);
      }
    }



    void analyze(const Event& e) {
      MSG_DEBUG("Num particles = " << apply<FinalState>(e, "FS").particles().size());

      const FastJets& jadejet = apply<FastJets>(e, "JadeJets");
      if (jadejet.clusterSeq()) {
        const double y_23 = jadejet.clusterSeq()->exclusive_ymerge_max(2);
        const double y_34 = jadejet.clusterSeq()->exclusive_ymerge_max(3);
        const double y_45 = jadejet.clusterSeq()->exclusive_ymerge_max(4);
        const double y_56 = jadejet.clusterSeq()->exclusive_ymerge_max(5);

        for (size_t i = 0; i < _h_R_Jade[0]->numBins(); ++i) {
          double ycut = _h_R_Jade[0]->bin(i).xMid();
          double width = _h_R_Jade[0]->bin(i).xWidth();
          if (y_23 < ycut) {
            _h_R_Jade[0]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Jade[1]->numBins(); ++i) {
          double ycut = _h_R_Jade[1]->bin(i).xMid();
          double width = _h_R_Jade[1]->bin(i).xWidth();
          if (y_34 < ycut && y_23 > ycut) {
            _h_R_Jade[1]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Jade[2]->numBins(); ++i) {
          double ycut = _h_R_Jade[2]->bin(i).xMid();
          double width = _h_R_Jade[2]->bin(i).xWidth();
          if (y_45 < ycut && y_34 > ycut) {
            _h_R_Jade[2]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Jade[3]->numBins(); ++i) {
          double ycut = _h_R_Jade[3]->bin(i).xMid();
          double width = _h_R_Jade[3]->bin(i).xWidth();
          if (y_56 < ycut && y_45 > ycut) {
            _h_R_Jade[3]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Jade[4]->numBins(); ++i) {
          double ycut = _h_R_Jade[4]->bin(i).xMid();
          double width = _h_R_Jade[4]->bin(i).xWidth();
          if (y_56 > ycut) {
            _h_R_Jade[4]->fillBin(i, width);
          }
        }
      }

      const FastJets& durjet = apply<FastJets>(e, "DurhamJets");
      if (durjet.clusterSeq()) {
        const double y_23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
        const double y_34 = durjet.clusterSeq()->exclusive_ymerge_max(3);
        const double y_45 = durjet.clusterSeq()->exclusive_ymerge_max(4);
        const double y_56 = durjet.clusterSeq()->exclusive_ymerge_max(5);

        _h_y_Durham[0]->fill(y_23);
        _h_y_Durham[1]->fill(y_34);
        _h_y_Durham[2]->fill(y_45);
        _h_y_Durham[3]->fill(y_56);

        for (size_t i = 0; i < _h_R_Durham[0]->numBins(); ++i) {
          double ycut = _h_R_Durham[0]->bin(i).xMid();
          double width = _h_R_Durham[0]->bin(i).xWidth();
          if (y_23 < ycut) {
            _h_R_Durham[0]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[1]->numBins(); ++i) {
          double ycut = _h_R_Durham[1]->bin(i).xMid();
          double width = _h_R_Durham[1]->bin(i).xWidth();
          if (y_34 < ycut && y_23 > ycut) {
            _h_R_Durham[1]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[2]->numBins(); ++i) {
          double ycut = _h_R_Durham[2]->bin(i).xMid();
          double width = _h_R_Durham[2]->bin(i).xWidth();
          if (y_45 < ycut && y_34 > ycut) {
            _h_R_Durham[2]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[3]->numBins(); ++i) {
          double ycut = _h_R_Durham[3]->bin(i).xMid();
          double width = _h_R_Durham[3]->bin(i).xWidth();
          if (y_56 < ycut && y_45 > ycut) {
            _h_R_Durham[3]->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[4]->numBins(); ++i) {
          double ycut = _h_R_Durham[4]->bin(i).xMid();
          double width = _h_R_Durham[4]->bin(i).xWidth();
          if (y_56 > ycut) {
            _h_R_Durham[4]->fillBin(i, width);
          }
        }
      }
    }



    /// Finalize
    void finalize() {
      for (size_t n = 0; n < 4; ++n) normalize(_h_y_Durham[n]);
      for (size_t n = 0; n < 5; ++n) scale(_h_R_Jade[n], 100/sumOfWeights());
      for (size_t n = 0; n < 5; ++n) scale(_h_R_Durham[n], 100/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_R_Jade[5];
    Histo1DPtr _h_R_Durham[5];
    Histo1DPtr _h_y_Durham[4];
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_OPAL_2000_S4300807);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// @brief CMS jet shape analysis
  /// @author Andreas Hinzmann
  class CMS_2012_I1111014 : public Analysis {
  public:

    /// Constructor
    CMS_2012_I1111014()
      : Analysis("CMS_2012_I1111014")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      const FinalState fs((Cuts::etaIn(-5.0, 5.0)));
      declare(fs, "FS");
      FastJets fj5(fs, FastJets::ANTIKT, 0.5);
      declare(fj5, "Jets5");
      FastJets fj7(fs, FastJets::ANTIKT, 0.7);
      declare(fj7, "Jets7");

      // Specify pT bins
      _ptedges = {{ 20.,25.,30.,40.,50.,60.,70.,80.,90.,100.,110.,125.,140.,160.,180.,200.,225.,250.,300.,400.,500.,600.,1000. }};
      _yedges  = {{ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 }};

      // Register a jet shape projection and histogram for each pT bin
      unsigned int histo_counter=1;
      for (size_t j = 0; j < 6; ++j) {
        for (size_t i = 0; i < 22; ++i) {
          if (i > 20 && j == 3) continue;
          if (i > 18 && j >= 4) continue;

          // Set up projections for each (pT,y) bin
          _jsnames_pT[i][j] = "JetShape" + to_str(i) + "_" + to_str(j);
          const JetShape jsp(fj7, 0.0, 0.7, 7, _ptedges[i], _ptedges[i+1], _yedges[j], _yedges[j+1], RAPIDITY);
          declare(jsp, _jsnames_pT[i][j]);

          // Book profile histograms for each (pT,y) bin
          book(_profhistRho_pT[i][j], histo_counter, 1, 1);
          histo_counter+=1;
        }
      }
      book(_profhistNch[0], 126, 1, 1);
      book(_profhistNch[1], 126, 1, 2);
      book(_profhistDr[0], 127, 1, 1);
      book(_profhistDr[1], 127, 1, 2);
      book(_profhistDeta, "TMP/Deta", refData(127,1,1));
      book(_profhistDphi, "TMP/Dphi", refData(127,1,1));
      book(_profhistAsym, "d128-x01-y01", true);

    }



    /// Do the analysis
    void analyze(const Event& evt) {

      // Get jets and require at least one to pass pT and y cuts
      Jets jets7 = apply<FastJets>(evt, "Jets7")
        .jetsByPt(Cuts::ptIn(_ptedges.front()*GeV, _ptedges.back()*GeV) && Cuts::absrap < 3.0);
      if(jets7.size()>2) jets7.resize(2); // Use only the first two jets
      MSG_DEBUG("Jet (R=0.7) multiplicity before cuts = " << jets7.size());
      if (jets7.size() == 0) {
        MSG_DEBUG("No jets (R=0.7) found in required pT and rapidity range");
        vetoEvent;
      }
      // Calculate and histogram jet shapes
      for (size_t jy = 0; jy < 6; ++jy) {
        for (size_t ipt = 0; ipt < 22; ++ipt) {
          if (ipt > 20 && jy == 3) continue;
          if (ipt > 18 && jy >= 4) continue;
          JetShape jsipt = apply<JetShape>(evt, _jsnames_pT[ipt][jy]);
          jsipt.calc(jets7);
          for (size_t ijet = 0; ijet < jsipt.numJets(); ++ijet) {
            for (size_t rbin = 0; rbin < jsipt.numBins(); ++rbin) {
              const double r_rho = jsipt.rBinMid(rbin);
              _profhistRho_pT[ipt][jy]->fill(r_rho, (1./0.1)*jsipt.diffJetShape(ijet, rbin));
            }
          }
        }
      }
      
      // Get jets and require at least one to pass pT and y cuts
      Jets jets5 = apply<FastJets>(evt, "Jets5")
        .jetsByPt(Cuts::ptIn(50*GeV, 1000*GeV) && Cuts::absrap < 2.0);
      // Calculate and histogram charged jet shapes
      for (const Jet& jet : jets5) {
        double ncharge = 0;
        double eta=0;
        double phi=0;
        double sumpt=0;
        for (const Particle& p : jet.particles()) {
          if ((p.pT() < 0.5) || (p.charge3()==0) || (abs(p.pid())==11) || (abs(p.pid())==13)) continue;
          ncharge++;
          sumpt+=p.pT();
          eta+=p.pT()*p.eta();
          phi+=p.pT()*mapAngleMPiToPi(p.phi()-jet.phi());
        }
        if (jet.absrap()<1.0) {
          _profhistNch[0]->fill(jet.pT(), ncharge );
        } else if (jet.absrap()<2.0) {
          _profhistNch[1]->fill(jet.pT(), ncharge );
        }
        if (sumpt==0) continue;
        eta/=sumpt;
        phi/=sumpt;
        double deta=0;
        double dphi=0;
        for (const Particle& p : jet.particles()) {
          if ((p.pT() < 0.5) || (p.charge3()==0) || (abs(p.pid())==11) || (abs(p.pid())==13)) continue;
          deta+=p.pT()*pow(p.eta()-eta,2);
          dphi+=p.pT()*pow(mapAngleMPiToPi(p.phi()-phi-jet.phi()),2);
        }
        deta/=sumpt;
        dphi/=sumpt;
        if ((dphi==0)||(deta==0)) continue;
        if (jet.absrap()<1.0) {
          _profhistDr[0]->fill(jet.pT(), deta+dphi );
          _profhistDeta->fill(jet.pT(), deta );
          _profhistDphi->fill(jet.pT(), dphi );
        } else if (jet.absrap()<2.0) {
          _profhistDr[1]->fill(jet.pT(), deta+dphi );
        }
      }
    }


    // Finalize
    void finalize() {
      for (unsigned int i = 0; i < _profhistAsym->numPoints(); ++i) {
        if((_profhistDeta->bin(i).effNumEntries()<2)||(_profhistDphi->bin(i).effNumEntries()<2)) continue;
        if((_profhistDeta->bin(i).mean()==0)||(_profhistDphi->bin(i).mean()==0)) continue;
        double mean_ratio=_profhistDeta->bin(i).mean() / _profhistDphi->bin(i).mean();
        double mean_error=mean_ratio*sqrt(pow(_profhistDeta->bin(i).stdErr()/_profhistDeta->bin(i).mean(),2)+pow(_profhistDphi->bin(i).stdErr()/_profhistDphi->bin(i).mean(),2));
        _profhistAsym->point(i).setY(mean_ratio,mean_error);
      }
    }

    //@}


  private:

    /// @name Analysis data
    //@{

    /// Jet \f$ p_\perp\f$ bins.
    vector<double> _ptedges; // This can't be a raw array if we want to initialise it non-painfully
    vector<double> _yedges;

    /// JetShape projection name for each \f$p_\perp\f$ bin.
    string _jsnames_pT[22][6];

    //@}

    /// @name Histograms
    //@{
    Profile1DPtr _profhistRho_pT[22][6];
    Profile1DPtr _profhistNch[2];
    Profile1DPtr _profhistDr[2];
    Profile1DPtr _profhistDeta;
    Profile1DPtr _profhistDphi;
    Scatter2DPtr _profhistAsym;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1111014);

}

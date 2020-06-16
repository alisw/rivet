// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief MC validation analysis for W polarisation
  class MC_WPOL : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_WPOL()
      : Analysis("MC_WPOL")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      WFinder wfinder(fs, Cuts::open(), PID::ELECTRON,
                      60.0*GeV, 100.0*GeV, 0.0*GeV, 0.0);
      declare(wfinder, "WFinder");
      Beam beams;
      declare(beams, "Beams");

      vector<string> tags{"_wplus", "_wminus"};
      _h_dists.resize(tags.size());
      _h_histos.resize(tags.size());
      for (size_t i=0; i<tags.size(); ++i) {
        _h_dists[i].resize(11,Profile1DPtr());
        double sqrts = sqrtS()>0. ? sqrtS() : 14000.;
        book(_h_dists[i][0] ,"A0"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][1] ,"A1"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][2] ,"A2"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][3] ,"A3"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][4] ,"A4"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][5] ,"A5"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][6] ,"A6"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][7] ,"A7"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][8] ,"fL"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][9] ,"fR"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        book(_h_dists[i][10] ,"f0"+tags[i],logspace(100, 1.0, 0.5*sqrts));
        _h_histos[i].resize(4,Histo1DPtr());
        book(_h_histos[i][0] ,"thetastar"+tags[i],100,-1.0,1.0);
        book(_h_histos[i][1] ,"phistar"+tags[i],90,0.0,360.0);
        book(_h_histos[i][2] ,"thetastar_ptw20"+tags[i],100,-1.0,1.0);
        book(_h_histos[i][3] ,"phistar_ptw20"+tags[i],90,0.0,360.0);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      const WFinder& wfinder = apply<WFinder>(event, "WFinder");
      if (wfinder.bosons().size() != 1) {
        vetoEvent;
      }
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();

      FourMomentum pb1(beams.second.momentum()), pb2(beams.first.momentum());
      Particle lepton = wfinder.constituentLeptons()[0];
      FourMomentum pl(lepton.momentum());
      size_t idx = (PID::charge3(lepton.pid())>0 ? 0 : 1);
      FourMomentum plnu(wfinder.bosons()[0].momentum());

      const LorentzTransform cms = LorentzTransform::mkFrameTransformFromBeta(plnu.betaVec());
      Matrix3 zrot(plnu.p3(), Vector3(0.0, 0.0, 1.0));
      pl=cms.transform(pl);
      pb1=cms.transform(pb1);
      pb2=cms.transform(pb2);
      Vector3 pl3=pl.p3();
      Vector3 pb13=pb1.p3();
      Vector3 pb23=pb2.p3();
      pl3=zrot*pl3;
      pb13=zrot*pb13;
      pb23=zrot*pb23;
      Vector3 xref(cos(pb13.theta())>cos(pb23.theta())?pb13:pb23);
      Matrix3 xrot(Vector3(xref.x(), xref.y(), 0.0), Vector3(1.0, 0.0, 0.0));
      pl3=xrot*pl3;

      double ptw(wfinder.bosons()[0].pT()/GeV);
      double thetas(pl3.theta()), phis(pl3.phi());
      double costhetas(cos(thetas)), sinthetas(sin(thetas));
      double cosphis(cos(phis)), sinphis(sin(phis));
      if (phis<0.0) phis+=2.0*M_PI;

      _h_histos[idx][0]->fill(costhetas,weight);
      _h_histos[idx][1]->fill(phis*180.0/M_PI,weight);
      if (ptw>20.0) {
        _h_histos[idx][2]->fill(costhetas,weight);
        _h_histos[idx][3]->fill(phis*180.0/M_PI,weight);
      }
      _h_dists[idx][0]->fill(ptw,10.0/3.0*(1.0-3.0*sqr(costhetas))+2.0/3.0,weight);
      _h_dists[idx][1]->fill(ptw,10.0*sinthetas*costhetas*cosphis,weight);
      _h_dists[idx][2]->fill(ptw,10.0*sqr(sinthetas)*(sqr(cosphis)-sqr(sinphis)),weight);
      _h_dists[idx][3]->fill(ptw,4.0*sinthetas*cosphis,weight);
      _h_dists[idx][4]->fill(ptw,4.0*costhetas,weight);
      _h_dists[idx][5]->fill(ptw,4.0*sinthetas*sinphis,weight);
      _h_dists[idx][6]->fill(ptw,10.0*costhetas*sinthetas*sinphis,weight);
      _h_dists[idx][7]->fill(ptw,10.0*sqr(sinthetas)*cosphis*sinphis,weight);
      _h_dists[idx][8]->fill(ptw,0.5*sqr(1.0-costhetas)-(1.0-2.0*sqr(costhetas)),weight);
      _h_dists[idx][9]->fill(ptw,0.5*sqr(1.0+costhetas)-(1.0-2.0*sqr(costhetas)),weight);
      _h_dists[idx][10]->fill(ptw,5.0*sqr(sinthetas)-3.0,weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for (size_t i=0; i<_h_histos.size(); ++i) {
        for (Histo1DPtr histo : _h_histos[i]) {
          scale(histo, crossSection()/picobarn/sumOfWeights());
        }
      }

    }

    //@}


  private:

    /// @name Histograms
    //@{

    vector<vector<Profile1DPtr> > _h_dists;
    vector<vector<Histo1DPtr> > _h_histos;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WPOL);

}

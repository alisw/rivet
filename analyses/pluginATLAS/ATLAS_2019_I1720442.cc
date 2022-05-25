// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// ATLAS 4-lepton lineshape at 13 TeV
  class ATLAS_2019_I1720442 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1720442);

    void init() {

      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState elecs(Cuts::abspid == PID::ELECTRON);
      PromptFinalState muons(Cuts::abspid == PID::MUON);

      // Selection
      Cut el_fid_sel = (Cuts::abseta < 2.47) && (Cuts::pT > 7*GeV);
      Cut mu_fid_sel = (Cuts::abseta < 2.7) && (Cuts::pT > 5*GeV);

      DressedLeptons dressed_elecs(photons, elecs, 0.005, el_fid_sel, false);
      declare(dressed_elecs, "elecs");

      DressedLeptons dressed_muons(photons, muons, 0.005, mu_fid_sel, false);
      declare(dressed_muons, "muons");

      // Book histos
      book(_h["m4l_inclusive"], 1,1,1);

      book(_h["m4l_ptslice1"], 2,1,1);
      book(_h["m4l_ptslice2"], 3,1,1);
      book(_h["m4l_ptslice3"], 4,1,1);
      book(_h["m4l_ptslice4"], 5,1,1);

      book(_h["m4l_rapidityslice1"], 6,1,1);
      book(_h["m4l_rapidityslice2"], 7,1,1);
      book(_h["m4l_rapidityslice3"], 8,1,1);
      book(_h["m4l_rapidityslice4"], 9,1,1);

      book(_h["m4l_4mu"], 12,1,1);
      book(_h["m4l_4e"], 13,1,1);
      book(_h["m4l_2e2mu"], 14,1,1);
    }


    /// @brief Generic dilepton candidate
    /// @todo Move into the Rivet core?
    struct Dilepton : public ParticlePair {
      Dilepton() { }
      Dilepton(ParticlePair _particlepair) : ParticlePair(_particlepair) {
        assert(first.abspid() == second.abspid());
      }
      FourMomentum mom() const { return first.momentum() + second.momentum(); }
      operator FourMomentum() const { return mom(); }
      static bool cmppT(const Dilepton& lx, const Dilepton& rx) { return lx.mom().pT() < rx.mom().pT(); }
      int flavour() const { return first.abspid(); }
      double pTl1() const { return first.pT(); }
      double pTl2() const { return second.pT(); }
    };


    struct Quadruplet {
      Quadruplet (Dilepton z1, Dilepton z2): _z1(z1), _z2(z2) { }
      enum class FlavCombi { mm=0, ee, me, em, undefined };
      FourMomentum mom() const { return _z1.mom() + _z2.mom(); }
      Dilepton getZ1() const { return _z1; }
      Dilepton getZ2() const { return _z2; }
      Dilepton _z1, _z2;
      FlavCombi type() const {
        if (     _z1.flavour() == 13 && _z2.flavour() == 13) { return FlavCombi::mm; }
        else if (_z1.flavour() == 11 && _z2.flavour() == 11) { return FlavCombi::ee; }
        else if (_z1.flavour() == 13 && _z2.flavour() == 11) { return FlavCombi::me; }
        else if (_z1.flavour() == 11 && _z2.flavour() == 13) { return FlavCombi::em; }
        else  return FlavCombi::undefined;
      }
    };


    vector<Quadruplet> getBestQuads(Particles& particles) {
      // H->ZZ->4l pairing
      // - Two same flavor opposite charged leptons
      // - Ambiguities in pairing are resolved by choosing the combination
      //     that results in the smaller value of |mll - mZ| for each pair successively
      vector<Quadruplet> quads {};

      size_t n_parts = particles.size();
      if (n_parts < 4)  return quads;

      // STEP 1: find SFOS pairs
      vector<Dilepton> SFOS;
      for (size_t i = 0; i < n_parts; ++i) {
        for (size_t j = 0; j < i; ++j) {
          if (particles[i].pid() == -particles[j].pid()) {
            // sort such that the negative lepton is listed first
            if (particles[i].pid() > 0)  SFOS.push_back(Dilepton(make_pair(particles[i], particles[j])));
            else                         SFOS.push_back(Dilepton(make_pair(particles[j], particles[i])));
          }
        }
      }
      if (SFOS.size() < 2)  return quads;

      // now we sort the SFOS pairs
      std::sort(SFOS.begin(), SFOS.end(), [](const Dilepton& p1, const Dilepton& p2) {
          return fabs(p1.mom().mass() - Z_mass) < fabs(p2.mom().mass() - Z_mass);
        });

      // Form all possible quadruplets passing the pT cuts
      for (size_t k = 0; k < SFOS.size(); ++k) {
        for (size_t l = k+1; l < SFOS.size(); ++l) {
          if (deltaR(SFOS[k].first.mom(),  SFOS[l].first.mom())  < 1e-13)  continue;
          if (deltaR(SFOS[k].first.mom(),  SFOS[l].second.mom()) < 1e-13)  continue;
          if (deltaR(SFOS[k].second.mom(), SFOS[l].first.mom())  < 1e-13)  continue;
          if (deltaR(SFOS[k].second.mom(), SFOS[l].second.mom()) < 1e-13)  continue;

          vector<double> lep_pt { SFOS[k].pTl1(), SFOS[k].pTl2(), SFOS[l].pTl1(), SFOS[l].pTl2() };
          std::sort(lep_pt.begin(), lep_pt.end(), std::greater<double>());
          if (!(lep_pt[0] > 20*GeV && lep_pt[1] > 15*GeV && lep_pt[2] > 10*GeV)) continue;
          quads.push_back( Quadruplet(SFOS[k], SFOS[l]) );
        }
      }
      return quads;
    }


    bool passMassCuts(const Quadruplet& theQuad){
      const vector<double> cuts_m34{ 5*GeV, 5*GeV, 12*GeV, 12*GeV, 50*GeV };
      const vector<double> cuts_m4l{ 0, 100*GeV, 110*GeV, 140*GeV, 190*GeV };

      double m4l = theQuad.mom().mass();
      double mZ1 = theQuad.getZ1().mom().mass();
      double mZ2 = theQuad.getZ2().mom().mass();

      // Invariant-mass requirements
      double cutval = cuts_m34.back();
      for (size_t k = 0; k < cuts_m34.size(); ++k) {
        if (cuts_m4l[k] > m4l) {
          cutval = cuts_m34[k-1] + (cuts_m34[k] - cuts_m34[k-1])/(cuts_m4l[k] - cuts_m4l[k-1]) * (m4l - cuts_m4l[k-1]);
          break;
        }
      }
      return inRange(mZ1, 50*GeV, 106*GeV) && inRange(mZ2, cutval, 115*GeV);
    }


    bool pass_dRll(const Quadruplet& theQuad) {
      const double dR_min_same = 0.1;
      const double dR_min_opp = 0.2;
      double dr_min_cross = dR_min_opp;
      if (theQuad.getZ1().flavour() == theQuad.getZ2().flavour()) {
        dr_min_cross = dR_min_same;
      }
      return !((deltaR(theQuad.getZ1().first,  theQuad.getZ1().second) < dR_min_same)  ||
               (deltaR(theQuad.getZ2().first,  theQuad.getZ2().second) < dR_min_same)  ||
               (deltaR(theQuad.getZ1().first,  theQuad.getZ2().first)  < dr_min_cross) ||
               (deltaR(theQuad.getZ1().first,  theQuad.getZ2().second) < dr_min_cross) ||
               (deltaR(theQuad.getZ1().second, theQuad.getZ2().first)  < dr_min_cross) ||
               (deltaR(theQuad.getZ1().second, theQuad.getZ2().second) < dr_min_cross));
    }


    bool pass_Jpsi(const Quadruplet& theQuad){
      Particles all_leps { theQuad.getZ1().first, theQuad.getZ1().second, theQuad.getZ2().first, theQuad.getZ2().second };
      for (const Particle& lep1 : all_leps) {
        for (const Particle& lep2 : all_leps) {
          if (lep1.pid() == -lep2.pid() && (lep1.mom() + lep2.mom()).mass() < 5*GeV) return false;
        }
      }
      return true;
    }


    // Handle 3 further CF stages - m12/34, dRmin, jpsi veto
    bool passSelection (const Quadruplet& theQuad){
      return passMassCuts(theQuad) && pass_dRll(theQuad) && pass_Jpsi(theQuad);
    }


    // Do the analysis
    void analyze(const Event& event) {

      //preselection of leptons for ZZ-> llll final state
      Particles dressed_leptons;
      for (auto lep : apply<DressedLeptons>(event, "muons").dressedLeptons()) { dressed_leptons.push_back(lep); }
      for (auto lep : apply<DressedLeptons>(event, "elecs").dressedLeptons()) { dressed_leptons.push_back(lep); }

      auto foundDressed = getBestQuads(dressed_leptons);
      // if we don't find any quad, we can stop here
      if (foundDressed.empty())  vetoEvent;

      bool pass = passSelection(foundDressed[0]);
      if (pass) {
        double m4l = foundDressed[0].mom().mass()/GeV;
        double pt4l = foundDressed[0].mom().pT()/GeV;
        double y4l = foundDressed[0].mom().absrap();
        Quadruplet::FlavCombi flavour = foundDressed[0].type();
        _h["m4l_inclusive"]->fill(m4l);
        if (     pt4l <  20.)  _h["m4l_ptslice1"]->fill(m4l);
        else if (pt4l <  50.)  _h["m4l_ptslice2"]->fill(m4l);
        else if (pt4l < 100.)  _h["m4l_ptslice3"]->fill(m4l);
        else if (pt4l < 600.)  _h["m4l_ptslice4"]->fill(m4l);

        if (     y4l < 0.4)  _h["m4l_rapidityslice1"]->fill(m4l);
        else if (y4l < 0.8)  _h["m4l_rapidityslice2"]->fill(m4l);
        else if (y4l < 1.2)  _h["m4l_rapidityslice3"]->fill(m4l);
        else if (y4l < 2.5)  _h["m4l_rapidityslice4"]->fill(m4l);

        if (     flavour == Quadruplet::FlavCombi::mm) _h["m4l_4mu"]->fill(m4l);
        else if (flavour == Quadruplet::FlavCombi::ee) _h["m4l_4e"]->fill(m4l);
        else if (flavour == Quadruplet::FlavCombi::me || flavour == Quadruplet::FlavCombi::em) {
          _h["m4l_2e2mu"]->fill(m4l);
        }
      }

    }


    /// Finalize
    void finalize() {
      const double sf = crossSection() / femtobarn / sumOfWeights();
      scale(_h, sf);
    }


  private:

    map<string, Histo1DPtr> _h;
    static constexpr double Z_mass = 91.1876;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1720442);

}

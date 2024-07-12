// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


    class MC_HFDECAYS : public Analysis {
    public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(MC_HFDECAYS);

      const string whoDis(const int pid) const {
        switch (pid) {
          case PID::B0:           return "B0";
          case PID::BPLUS:        return "BPLUS";
          case PID::B0S:          return "B0S";
          case PID::LAMBDAB:      return "LAMBDAB";
          case PID::D0:           return "D0";
          case PID::DPLUS:        return "DPLUS";
          case PID::DSPLUS:       return "DSPLUS";
          case PID::LAMBDACPLUS:  return "LAMBDACPLUS";
          default:                return "";
        }
      }

      void addBinned(const string &name, const YODA::HistoBin1D &ptbin, const int nbins, const double lo, const double hi) {
        string suff = "_" + to_string(int(ptbin.xMin())) + "_" + to_string(int(ptbin.xMax()));
        { Histo1DPtr tmp; _b[name].add(ptbin.xMin(), ptbin.xMax(), book(tmp, name+suff, nbins, lo, hi)); }
      }

      /// Book histograms and initialise projections before the run
      void init() {

        declare(HeavyHadrons(),"HA");

        FastJets jetpro(FinalState(), FastJets::ANTIKT, 0.4, JetAlg::Muons::DECAY, JetAlg::Invisibles::DECAY);
        declare(jetpro, "Jets");

        // bar charts
        book(_h["bar_b_jet_width"], "width_B_jet", 7, 0., 0.3);
        book(_h["bar_c_jet_width"], "width_C_jet", 7, 0., 0.3);
        // profiles
        book(_p["b_jet_rho"], "avg_rho_B_jet", 10, 0., 0.4);
        book(_p["c_jet_rho"], "avg_rho_C_jet", 10, 0., 0.4);
        book(_p["b_jet_psi"], "avg_psi_B_jet", 10, 0., 0.4);
        book(_p["c_jet_psi"], "avg_psi_C_jet", 10, 0., 0.4);
        // histograms
        book(_h["b_frac"], "prod_frac_B", 11, 0.5,11.5);
        book(_h["c_frac"], "prod_frac_C",  7, 0.5, 7.5);
        book(_h["b_jet_frag"], "frag_B_jet", 50, 0., 1.1);
        book(_h["c_jet_frag"], "frag_C_jet", 50, 0., 1.1);
        book(_h["b_jet_pT"], "pT_B_jet", 25, 25., 1025.);
        book(_h["c_jet_pT"], "pT_C_jet", 25, 25., 1025.);
        book(_h["b_jet_pThad"], "pThad_B_jet", 10, 0., 200.);
        book(_h["c_jet_pThad"], "pThad_C_jet", 10, 0., 200.);
        book(_h["b_jet_high_pT"], "high_pT_frag_B_jet", 46, 0., 1.);
        book(_h["c_jet_high_pT"], "high_pT_frag_C_jet", 46, 0., 1.);
        book(_h["b_jet_high_pT_1H"], "high_pT_frag_B_jet_1H", 46, 0., 1.);
        book(_h["c_jet_high_pT_1H"], "high_pT_frag_C_jet_1H", 46, 0., 1.);


        book(_h["ch_B0"         ], "B0_charged_mult",          30, 0.5, 30.5);
        book(_h["ch_BPLUS"      ], "BPLUS_charged_mult",       30, 0.5, 30.5);
        book(_h["ch_B0S"        ], "B0S_charged_mult",         30, 0.5, 30.5);
        book(_h["ch_LAMBDAB"    ], "LAMBDAB_charged_mult",     30, 0.5, 30.5);
        book(_h["ch_D0"         ], "D0_charged_mult",          16, 0.5, 16.5);
        book(_h["ch_DPLUS"      ], "DPLUS_charged_multh",      16, 0.5, 16.5);
        book(_h["ch_DSPLUS"     ], "DSPLUS_charged_mult",      16, 0.5, 16.5);
        book(_h["ch_LAMBDACPLUS"], "LAMBDACPLUS_charged_mult", 16, 0.5, 16.5);

        book(_h["st_B0"         ], "B0_stable_mult",          40, 0.5, 40.5);
        book(_h["st_BPLUS"      ], "BPLUS_stable_mult",       40, 0.5, 40.5);
        book(_h["st_B0S"        ], "B0S_stable_mult",         40, 0.5, 40.5);
        book(_h["st_LAMBDAB"    ], "LAMBDAB_stable_mult",     40, 0.5, 40.5);
        book(_h["st_D0"         ], "D0_stable_mult",          20, 0.5, 20.5);
        book(_h["st_DPLUS"      ], "DPLUS_stable_mult",       20, 0.5, 20.5);
        book(_h["st_DSPLUS"     ], "DSPLUS_stable_mult",      20, 0.5, 20.5);
        book(_h["st_LAMBDACPLUS"], "LAMBDACPLUS_stable_mult", 20, 0.5, 20.5);

        book(_h["pt_B0"         ], "B0_pT",          10, 25., 425.);
        book(_h["pt_BPLUS"      ], "BPLUS_pT",       10, 25., 425.);
        book(_h["pt_B0S"        ], "B0S_pT",         10, 25., 425.);
        book(_h["pt_LAMBDAB"    ], "LAMBDAB_pT",     10, 25., 425.);
        book(_h["pt_D0"         ], "D0_pT",          10, 25., 425.);
        book(_h["pt_DPLUS"      ], "DPLUS_pT",       10, 25., 425.);
        book(_h["pt_DSPLUS"     ], "DSPLUS_pT",      10, 25., 425.);
        book(_h["pt_LAMBDACPLUS"], "LAMBDACPLUS_pT", 10, 25., 425.);

        book(_h["b_jet_ch_mult"], "charged_mult_B_jets", 40, 0.5, 40.5);
        book(_h["c_jet_ch_mult"], "charged_mult_C_jets", 40, 0.5, 40.5);

        book(_h["b_jet_l_pTrel"], "lepton_pTrel_B_jets", 8, 0., 15.);
        book(_h["c_jet_l_pTrel"], "lepton_pTrel_C_jets", 8, 0., 10.);

        book(_h["b_jet_l_pT"], "lepton_pT_B_jets", 10, 0., 100.);
        book(_h["c_jet_l_pT"], "lepton_pT_C_jets", 10, 0., 100.);

        // double-differentials
        ptaxis = YODA::Histo1DAxis({25, 30, 50, 70, 100, 150, 300, 500, 1000});
        for (size_t i = 0; i < ptaxis.numBins(); ++i) {
          string suff = to_string(int(ptaxis.bin(i).xMin())) + "_" + to_string(int(ptaxis.bin(i).xMax()));
          book(_p["avg_B_jet_rho_"+to_string(i)], "avg_B_jet_rho_"+suff, 10, 0., 0.4);
          book(_p["avg_C_jet_rho_"+to_string(i)], "avg_C_jet_rho_"+suff, 10, 0., 0.4);

          addBinned("avg_B_jet_ch_mult", ptaxis.bin(i), 40, 0.5, 40.5);
          addBinned("avg_C_jet_ch_mult", ptaxis.bin(i), 40, 0.5, 40.5);
          if (i == 0) {
            addBinned("avg_B_jet_l_pTrel", ptaxis.bin(i), 2, 0., 4.);
            addBinned("avg_C_jet_l_pTrel", ptaxis.bin(i), 2, 0., 3.);
          }
          else if (i <= 2) {
            addBinned("avg_B_jet_l_pTrel", ptaxis.bin(i), 4, 0., 8.);
            addBinned("avg_C_jet_l_pTrel", ptaxis.bin(i), 5, 0., 6.);
          }
          else {
            addBinned("avg_B_jet_l_pTrel", ptaxis.bin(i), 8, 0., 15.);
            addBinned("avg_C_jet_l_pTrel", ptaxis.bin(i), 8, 0., 10.);
          }
        }
      }

      double p_annulus(const Jet &jet, const double a, const double b) const {
        // calculate the total momentum inside an annulus with a <= R < b
        return sum(filter_select(jet.particles(), [&](const Particle &p) {
          const double dr = deltaR(p, jet);
          return (dr < b && dr >= a);
        }), Kin::pT, 0.)/GeV;
      }

      void count_mult(const Particle& p) {
        unsigned int nst = p.stableDescendants().size() + 0.5;
        unsigned int nch = p.stableDescendants(Cuts::charge != 0).size() + 0.5;
        _h["st_"+whoDis(p.pid())]->fill(nst);
        _h["ch_"+whoDis(p.pid())]->fill(nch);
      }

      double pTrel (const Jet& jet, const Particle& p) const {
        return (p.p3().cross(jet.p3())).mod()/(jet.p3().mod());
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        const HeavyHadrons &ha = apply<HeavyHadrons>(event,"HA");

        if (ha.bHadrons().empty() && ha.cHadrons().empty()) vetoEvent;

        for (const Particle &p : ha.bHadrons()) {
          const string name = "pt_" + whoDis(p.pid());
          if (p.pid() == PID::B0) {
            count_mult(p);
            _h["b_frac"]->fill(1);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::BPLUS) {
            count_mult(p);
            _h["b_frac"]->fill(2);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::B0S) {
            count_mult(p);
            _h["b_frac"]->fill(3);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::BCPLUS)   _h["b_frac"]->fill(4);
          else if (p.pid() == PID::LAMBDAB) {
            count_mult(p);
            _h["b_frac"]->fill(5);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::XIBMINUS)     _h["b_frac"]->fill(6);
          else if (p.pid() == PID::XI0B)         _h["b_frac"]->fill(7);
          else if (p.pid() == PID::OMEGABMINUS)  _h["b_frac"]->fill(8);
          else if (p.pid() == PID::SIGMABMINUS)  _h["b_frac"]->fill(9);
          else if (p.pid() == PID::SIGMAB)       _h["b_frac"]->fill(10);
          else if (p.pid() == PID::SIGMABPLUS)   _h["b_frac"]->fill(11);
        }

        for (const Particle &p : ha.cHadrons()) {
          const string name = "pt_" + whoDis(p.pid());
          if (p.pid() == PID::DPLUS) {
            count_mult(p);
            _h["c_frac"]->fill(1);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::D0) {
            count_mult(p);
            _h["c_frac"]->fill(2);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::DSPLUS) {
            count_mult(p);
            _h["c_frac"]->fill(3);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::LAMBDACPLUS) {
            count_mult(p);
            _h["c_frac"]->fill(4);
            _h[name]->fill(p.pT()/GeV);
          }
          else if (p.pid() == PID::XICPLUS)  _h["c_frac"]->fill(5);
          else if (p.pid() == PID::XI0C)     _h["c_frac"]->fill(6);
          else if (p.pid() == PID::OMEGA0C)  _h["c_frac"]->fill(7);
        }

        Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25.*GeV && Cuts::absrap < 2.5);
        if (jets.empty()) vetoEvent;

        auto r_bins = _p["b_jet_rho"]->xEdges();
        const double dr = r_bins[1]-r_bins[0]; //dr is equal for all bins
        for (const Jet& thisJet : jets) {
          double p_0_R = p_annulus(thisJet, 0., 0.4);
          if (fuzzyEquals(p_0_R, 0., 1e-5))  continue;
          Particles bjets = thisJet.bTags(Cuts::pT > 5*GeV);
          ifilter_select(bjets, deltaRLess(thisJet, 0.3));
          for (const Particle &thisB : bjets) {
            _h["b_jet_pThad"]->fill(thisB.pT()/GeV);
            double z = thisJet.p3().dot(thisB.p3())/thisJet.p2();
            _h["b_jet_frag"]->fill(z);
            if (inRange(thisJet.pT(), 500*GeV, 1000*GeV)) {
              _h["b_jet_high_pT"]->fill(z);
              if (bjets.size() == 1)  _h["b_jet_high_pT_1H"]->fill(z);
            }
            for (size_t i = 0; i < r_bins.size()-1; ++i) {
              double r = 0.5*(r_bins[i]+r_bins[i+1]);
              double rho = p_annulus(thisJet, r_bins[i], r_bins[i+1])/p_0_R/dr;
              double psi = p_annulus(thisJet, 0, r_bins[i+1])/p_0_R;
              _p["b_jet_rho"]->fill(r, rho);
              _p["b_jet_psi"]->fill(r, psi);
              size_t ptb = ptaxis.binIndexAt(thisJet.pT()/GeV); // index in {0, nBins - 1}
              if (ptb < ptaxis.numBins()) _p["avg_B_jet_rho_"+to_string(ptb)]->fill(r, rho);
            }
          }

          Particles cjets = thisJet.cTags(Cuts::pT > 5*GeV);
          ifilter_select(cjets, deltaRLess(thisJet, 0.3));
          if (bjets.empty()) {
            for (const Particle& thisC : cjets) {
              _h["c_jet_pThad"]->fill(thisC.pT()/GeV);
              double z = thisJet.p3().dot(thisC.p3())/thisJet.p2();
              _h["c_jet_frag"]->fill(z);
              if (inRange(thisJet.pT(), 500*GeV, 1000*GeV)) {
                _h["c_jet_high_pT"]->fill(z);
                if (cjets.size() == 1) {
                  _h["c_jet_high_pT_1H"]->fill(z);
                }
              }
              for (size_t i = 0; i < r_bins.size()-1; ++i) {
                double r = 0.5*(r_bins[i]+r_bins[i+1]);
                double rho = p_annulus(thisJet, r_bins[i], r_bins[i+1])/p_0_R/dr;
                double psi = p_annulus(thisJet, 0, r_bins[i+1])/p_0_R;
                _p["c_jet_rho"]->fill(r, rho);
                _p["c_jet_psi"]->fill(r, psi);
                size_t ptb = ptaxis.binIndexAt(thisJet.pT()/GeV); // index in {0, nBins - 1}
                if (ptb < ptaxis.numBins()) _p["avg_C_jet_rho_"+to_string(ptb)]->fill(r, rho);
              }
            }
          }

          if (bjets.size()) {
            double W_num = 0., W_den = 0.;
            long N_charged = 0;
            for (const Particle& pp : thisJet.particles()) {
              if(pp.isStable()) {
                W_num += pp.pT()*deltaR(thisJet,pp);
                W_den += pp.pT();
                if (pp.isCharged())  ++N_charged;
              }
              if(pp.isLepton()) {
                _h["b_jet_l_pT"]->fill(pp.pT()/GeV);
                _h["b_jet_l_pTrel"]->fill(pTrel(thisJet,pp)/GeV);
                _b["avg_B_jet_l_pTrel"].fill(thisJet.pT()/GeV, pTrel(thisJet,pp)/GeV);
              }
            }
            if (W_den)  _h["bar_b_jet_width"]->fill(W_num/W_den);
            _h["b_jet_pT"]->fill(thisJet.pT()/GeV);
            if (N_charged) {
              _h["b_jet_ch_mult"]->fill(N_charged);
              _b["avg_B_jet_ch_mult"].fill(thisJet.pT()/GeV, N_charged);
            }
          }
          else if(cjets.size()) {
            double W_num = 0., W_den = 0.;
            long N_charged = 0;
            for (const Particle& pp : thisJet.particles()) {
              if(pp.isStable()) {
                W_num += pp.pT()*deltaR(thisJet,pp);
                W_den += pp.pT();
                if(pp.isCharged()) ++N_charged;
              }
              if(pp.isLepton()) {
                _h["c_jet_l_pT"]->fill(pp.pT()/GeV);
                _h["c_jet_l_pTrel"]->fill(pTrel(thisJet,pp)/GeV);
                _b["avg_C_jet_l_pTrel"].fill(thisJet.pT()/GeV, pTrel(thisJet,pp)/GeV);
              }
            }
            if (W_den)  _h["bar_c_jet_width"]->fill(W_num/W_den);
            _h["c_jet_pT"]->fill(thisJet.pT()/GeV);
            if (N_charged) {
              _h["c_jet_ch_mult"]->fill(N_charged);
              _b["avg_C_jet_ch_mult"].fill(thisJet.pT()/GeV, N_charged);
            }
          }
        }
      }

      /// Normalise histograms etc., after the run
      void finalize() {
        for (const auto &hit : _h) {
          double sf = 1.0;
          if (hit.first.find("bar_") != string::npos)  sf = (hit.second->xMax()-hit.second->xMin())/hit.second->numBins();
          normalize(hit.second, sf);
        }
        for (const auto &bit : _b) {
          for (const auto& hist : bit.second.histos()) { normalize(hist); }
        }
      }

    private:
      map<string, Histo1DPtr> _h;
      map<string, Profile1DPtr> _p;
      map<string, BinnedHistogram> _b;
      YODA::Histo1DAxis ptaxis;

  };

  RIVET_DECLARE_PLUGIN(MC_HFDECAYS);
}

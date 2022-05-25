// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {

/// @brief Measurements of differential production cross sections for a Z boson in association with jets in pp collisions at 8 TeV

  class CMS_2017_I1497519 : public Analysis {
  private:
    enum histIds {
      //Normalized differential cross sections
      kYZ, kYZmidPt, kYZhighPt, //Figs. 9a, 9b, 9c
      kYdiff, kYdiffMidPt, kYdiffHighPt, //Figs. 10a, 10b, 10c
      kYdiff2jZj1, kYdiff2jZj2, kYdiff2jZdijet, kYdiff2jJJ, //Figs. 11a, 11b, 11c, 12a
      kYsum, kYsumMidPt, kYsumHighPt, //Figs. 13a, 13b, 13c
      kYsum2jZj1, kYsum2jZj2, kYsum2jZdijet, kYsum2jJJ, //Figs. 14a, 14b, 14c, 12b

      //Absolute differential cross sections. Figure numbers refer to doi:10.1007/JHEP04(2017)022
      kNjets_exc, kNjets_inc, //Figs. 2a, 2b
      kPtj1, kPtj2, kPtj3, kPtj4, kPtj5, //Figs. 3a, 3b, 4a, 4b, 5
      kYj1, kYj2, kYj3, kYj4, kYj5, //Figs. 6a, 6b, 7a, 7b, 8
      kHt1, kHt2, kHt3, kHt4, kHt5, //Figs. 15a, 15b, 16a, 16b, 17
      kDphiZj1, kDphiZj1MidPt, kDphiZj1HighPt, // Figs. 18a, 20a, 22a
      kDphi2jZj1, kDphi2jZj1MidPt,  kDphi2jZj1HighPt, //Figs. 18b, 20b, 22b
      kDphi3jZj1, kDphi3jZj1MidPt, kDphi3jZj1HighPt,// Figs. 18c, 20c, 22c
      kDphi3jZj2, kDphi3jZj2MidPt, kDphi3jZj2HighPt, // Figs. 19a, 21a, 23a
      kDphi3jZj3, kDphi3jZj3MidPt, kDphi3jZj3HighPt, // Figs. 19b, 21b, 23b
      kDphi3jZj1HighHT, kDphi3jZj2HighHT, kDphi3jZj3HighHT, // Figs. 24a, 24b, 24c
      kDphi3jJ1j2, kDphi3jJ1j2MidPt, kDphi3jJ1j2HighPt, //Figs. 25a, 26a, 27a
      kDphi3jJ1j3, kDphi3jJ1j3MidPt, kDphi3jJ1j3HighPt, //Figs. 25b, 26b, 27b
      kDphi3jJ2j3, kDphi3jJ2j3MidPt, kDphi3jJ2j3HighPt, //Figs. 25c, 26c, 27c
      kDijetMass, //Fig. 28
      kPtjYjBin0, kPtjYjBin1, kPtjYjBin2, kPtjYjBin3, kPtjYjBin4, kPtjYjBin5, kPtjYjBin6, //Fig. 29
      kYjYZBin0, kYjYZBin1, kYjYZBin2, kYjYZBin3, //Fig. 33
      kPtjYjYZBin0, kPtjYjYZBin1, kPtjYjYZBin2, //Fig. 37 top left frame (Z and j on same side)
      kPtjYjYZBin3, kPtjYjYZBin4, kPtjYjYZBin5, //Fig. 37 bottom left frame  (same side)
      kPtjYjYZBin6, kPtjYjYZBin7, kPtjYjYZBin8, //Fig. 37 top right frame (opposite sides)
      kPtjYjYZBin9, kPtjYjYZBin10, kPtjYjYZBin11, //Fig. 37 bottom right frame (opposite sides)
      nHistos
    };

    unsigned nNormalized = kYsum2jJJ + 1;

    template<typename T>
    inline double ydiff(const ParticleBase& p1, const T& p2){ return 0.5*fabs(p1.rapidity()-p2.rapidity()); }

    template<typename T>
    inline double ysum(const ParticleBase&  p1, const T& p2){ return 0.5*fabs(p1.rapidity()+p2.rapidity()); }

    /// Fill histograms _h[ih+0..2] with x with different lower thresholds on tocut:
    /// no threshold, threshold cut2, threshold cut2
    void fill3cuts(int ih, double tocut, double cut2, double cut3, double x){
      _h[ih]->fill(x);
      if (tocut > cut2) _h[ih+1]->fill(x);
      if (tocut > cut3) _h[ih+2]->fill(x);
    }

    /// Fills a set of N 1-D histograms _h[ih + 0...N-1] that represents a 2-D distribution
    /// @param bins: boundaries of the y bins covered by each 1-D histogram
    void fill2D(int ih, std::vector<double>& ybins, double x, double y, double w){
      int iybin = -1;
      double lowEdge;
      for(auto highEdge: ybins){
        if(y < highEdge){
          if (iybin >= 0) _h[ih + iybin]->fill(x, w / (highEdge - lowEdge));
          break;
        }
        lowEdge = highEdge;
        ++iybin;
      }
    }

    void fill3D(int ih, std::vector<double>& ybins, std::vector<double> zbins, double x, double y, double z, double w){
      int izbin = -1;
      double lowEdge = 0;
      int nybins = ybins.size() - 1;
      for(auto highEdge: zbins){
        if(z < highEdge){
          if (izbin >= 0) fill2D(ih + izbin*nybins, ybins, x, y, w / (highEdge - lowEdge));
          break;
        }
        lowEdge = highEdge;
        ++izbin;
      }
    }

  public:

    /// Constructor
    CMS_2017_I1497519()
      : Analysis("CMS_2017_I1497519")
    {
      _histListInPaperOrder =  { /*number NN in comment is the id from the histogram name dNN-x01-y01*/
        /*1*/kNjets_exc, /*2*/ kNjets_inc, //Figs. 2a, 2b
        /*3*/kPtj1, /*4*/kPtj2, /*5*/kPtj3, /*6*/kPtj4, /*7*/kPtj5, //Figs. 3a, 3b, 4a, 4b, 5
        /*8*/kYj1, /*9*/kYj2, /*10*/kYj3, /*11*/kYj4, /*12*/kYj5, //Figs. 6a, 6b, 7a, 7b, 8
        /*13*/kYZ, /*14*/kYZmidPt, /*15*/kYZhighPt, //Figs. 9a, 9b, 9c
        /*16*/kYdiff, /*17*/kYdiffMidPt, /*18*/kYdiffHighPt, //Figs. 10a, 10b, 10c
        /*19*/kYdiff2jZj1, /*20*/kYdiff2jZj2, /*21*/kYdiff2jZdijet, /*22*/kYdiff2jJJ, //Figs. 11a, 11b, 11c, 12a
        /*23*/kYsum2jJJ, //12b
        /*24*/kYsum, /*25*/kYsumMidPt, /*26*/kYsumHighPt, //Figs. 13a, 13b, 13c
        /*27*/kYsum2jZj1, /*28*/kYsum2jZj2, /*29*/kYsum2jZdijet, //Figs. 14a, 14b, 14c
        /*30*/kHt1, /*31*/kHt2, /*32*/kHt3, /*33*/kHt4, /*34*/kHt5, //Figs. 15a, 15b, 16a, 16b, 17
        /*35*/kDphiZj1, /*36*/kDphi2jZj1, /*37*/kDphi3jZj1,  // Figs. 18a, 18b, 18c
        /*38*/kDphi3jZj2, /*39*/kDphi3jZj3, // Figs. 19a, 19b
        /*40*/kDphiZj1MidPt, /*41*/kDphi2jZj1MidPt, /*42*/kDphi3jZj1MidPt, // Figs. 20a, 20b, 20c
        /*43*/kDphi3jZj2MidPt, /*44*/kDphi3jZj3MidPt, // Figs. 21a, 21b
        /*45*/kDphiZj1HighPt, /*46*/kDphi2jZj1HighPt, /*47*/kDphi3jZj1HighPt,// Figs. 22a, 22b, 22c
        /*48*/kDphi3jZj2HighPt, /*49*/kDphi3jZj3HighPt, // Figs. 23a, 23b
        /*50*/kDphi3jZj1HighHT, /*51*/kDphi3jZj2HighHT, /*52*/kDphi3jZj3HighHT, // Figs. 24a, 24b, 24c
        /*53*/kDphi3jJ1j2, /*54*/kDphi3jJ1j3, /*55*/kDphi3jJ2j3, //Figs. 25a, 25b, 25c
        /*56*/kDphi3jJ1j2MidPt, /*57*/kDphi3jJ1j3MidPt, /*58*/kDphi3jJ2j3MidPt, //Figs. 26a, 26b, 26c
        /*59*/kDphi3jJ1j2HighPt, /*60*/kDphi3jJ1j3HighPt, /*61*/kDphi3jJ2j3HighPt, //Figs. 27a, 27b, 27c
        /*62*/kDijetMass, //Fig. 28
        /*63*/kPtjYjBin0, /*64*/kPtjYjBin1, /*65*/kPtjYjBin2, /*66*/kPtjYjBin3, /*67*/kPtjYjBin4, /*68*/kPtjYjBin5, /*69*/kPtjYjBin6, //Fig. 29
        /*70*/kYjYZBin0, /*71*/kYjYZBin1, /*72*/kYjYZBin2, /*73*/kYjYZBin3, //Fig. 33
        /*74*/kPtjYjYZBin0, /*75*/kPtjYjYZBin1, /*76*/kPtjYjYZBin2, //Fig. 37 top left frame (Z and j on same side)
        /*77*/kPtjYjYZBin6, /*78*/kPtjYjYZBin7, /*79*/kPtjYjYZBin8, //Fig. 37 top right frame (opposite sides)
        /*80*/kPtjYjYZBin3, /*81*/kPtjYjYZBin4, /*82*/kPtjYjYZBin5, //Fig. 37 bottom left frame  (same side)
        /*83*/kPtjYjYZBin9, /*84*/kPtjYjYZBin10, /*85*/kPtjYjYZBin11, //Fig. 37 bottom right frame (opposite sides)
      };
    }


    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      // default to combined.
      _mode = 2;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;
      if ( getOption("LMODE") == "EMU" ) _mode = 2;

      FinalState fs;

      PromptFinalState bareMuons(Cuts::abspid == PID::MUON);
      declare(DressedLeptons(fs, bareMuons, /*dRmax = */0.1,
                             Cuts::pT > 20*GeV && Cuts::abseta < 2.4,
                             /*useDecayPhotons = */ true),
              "muons");

      PromptFinalState bareElectrons(Cuts::abspid == PID::ELECTRON);
      declare(DressedLeptons(fs, bareElectrons, /*dRmax =*/ 0.1,
                             Cuts::pT > 20*GeV && Cuts::abseta < 2.4,
                             /*useDecayPhotons = */ true),
              "electrons");

      FastJets jets(fs, FastJets::ANTIKT, 0.5);
      declare(jets, "jets");

      _h = std::vector<Histo1DPtr>(nHistos);
      for(int ih = 0; ih < nHistos; ++ih){
        book(_h[_histListInPaperOrder[ih]], ih + 1, 1, 1);
      }

      _ptjYjBins = {0, 0.5, 1., 1.5, 2., 2.5, 3.2, 4.7};
      _yJyZbins = {0, 0.5, 1., 1.5, 2.5};
      _ptJyJyZbinsYj = {0., 1.5, 2.5, 4.7};
      _ptJyJyZbinsYZ = {0., 1., 2.5};
    }

    /// Z boson finder.
    /// Note: we don't use the standard ZFinder class in order to stick to
    /// the definition of the publication that is simpler than the ZFinder
    /// algorithm
    /// @param leptons pt-ordered of electron or muon collection to use to build
    /// the Z boson
    std::unique_ptr<Particle> zfinder(const std::vector<DressedLepton>& leptons){
      if(leptons.size() < 2) return 0;
      if(leptons[0].charge()*leptons[1].charge() > 0) return 0;
      std::unique_ptr<Particle> cand(new Particle(PID::ZBOSON, leptons[0].mom()
                                                  + leptons[1].mom()));
      if (cand->mass() < 71.*GeV || cand->mass() > 111.*GeV) return 0;
      return cand;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      std::vector<DressedLepton> muons = apply<DressedLeptons>(event, "muons").dressedLeptons();
      std::vector<DressedLepton> electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();

      //Look for Z->ee
      std::unique_ptr<Particle> z = zfinder(electrons);

      const std::vector<DressedLepton>* dressedLeptons = 0;

      //Look for Z->ee
      if(z.get() != nullptr && _mode != 1) {
        dressedLeptons = &electrons;
      } else{ //look for Z->mumu
        z = zfinder(muons);
        if(z.get() != nullptr && _mode != 0){
          dressedLeptons = &muons;
        } else{ //no Z boson found
          vetoEvent;
        }
      }

      // Cluster jets
      const FastJets& fj = apply<FastJets>(event, "jets");
      const Jets& jets = fj.jetsByPt(Cuts::absrap < 4.7 && Cuts::pT > 30*GeV);

      // Remove jets overlapping with any of the two selected leptons
      Jets goodjets47 = filter_discard(jets, [dressedLeptons](const ParticleBase& j){
          return deltaR(j, (*dressedLeptons)[0]) < 0.5
          ||  deltaR(j, (*dressedLeptons)[1]) < 0.5;
        });

      // Jets in the CMS tracker acceptance
      Jets goodjets24 = filter_select(goodjets47, [](const ParticleBase& j){
          return j.absrapidity() < 2.4;
        });

      goodjets24 = sortByPt(goodjets24);

      // Compute jet pt scalar sum, H_T:
      double ht = sum(goodjets24, Kin::pT, 0.);

      _h[kNjets_exc]->fill(goodjets24.size());

      // Fill jet number integral histograms
      /// @todo Could be better computed by toIntegral transform on exclusive histo
      for (size_t iJet = 0; iJet <= goodjets24.size(); iJet++ )
        _h[kNjets_inc]->fill(iJet);

      int offset = 0;
      for(const auto& j: goodjets24){
        int nj = 1 + offset;
        if(nj > 5) break;

        _h[kPtj1 + offset]->fill(j.pT() / GeV);
        _h[kYj1 + offset]->fill(j.absrapidity());
        _h[kHt1 + offset]->fill(ht);
        ++offset;
      }

      ///////////////////////////////
      /// Nj >=1  in |y| < 4.7
      if (goodjets47.size() < 1) return;

      const Jet& j1_47 = goodjets47[0];

      //note: 0.5 factors in the four following fill statements is required
      //because the cross section is differiated in y(j1) while the binning
      //is in abs(y(j1)): bin filled twice for the same y(j1) value.
      //An extra factor 0.5 is included for differential cross including jet and
      //Z boson rapidities to matches the definition used in the publication.
      fill2D(kPtjYjBin0, _ptjYjBins, j1_47.pt(), j1_47.absrapidity(), 0.5);
      fill2D(kYjYZBin0, _yJyZbins, j1_47.absrapidity()*sign(z->rapidity()*j1_47.rapidity()), z->absrapidity(), 0.25);

      if(j1_47.rapidity()*z->rapidity() >= 0){
        fill3D(kPtjYjYZBin0, _ptJyJyZbinsYj, _ptJyJyZbinsYZ, j1_47.pt(), j1_47.absrapidity(), z->absrapidity(), 0.25);
      } else{
        fill3D(kPtjYjYZBin6, _ptJyJyZbinsYj, _ptJyJyZbinsYZ, j1_47.pt(), j1_47.absrapidity(), z->absrapidity(), 0.25);
      }
      ///////////////////////////////

      ////////////////////////////////
      /// Nj >= 1 in |y| < 2.4
      if (goodjets24.size() < 1) return;

      const Jet& j1 = goodjets24[0];

      fill3cuts(kYZ, z->pt(), 150*GeV, 300*GeV, z->absrapidity());

      fill3cuts(kYdiff, z->pt(), 150*GeV, 300*GeV, ydiff(*z, j1));
      fill3cuts(kYsum, z->pt(), 150*GeV, 300*GeV, ysum(*z, j1));

      fill3cuts(kDphiZj1, z->pt(), 150*GeV, 300*GeV, deltaPhi(*z, j1));
      //////////////////////////////

      ////////////////////////////////
      /// Nj >= 2 in |y| < 2.4
      if (goodjets24.size() < 2) return;

      const Jet& j2 = goodjets24[1];
      _h[kYdiff2jZj1]->fill(ydiff(*z, j1));
      _h[kYdiff2jZj2]->fill(ydiff(*z, j2));
      _h[kYdiff2jZdijet]->fill(ydiff(*z, j1.mom() + j2.mom()));
      _h[kYdiff2jJJ]->fill(ydiff(j1, j2));

      _h[kYsum2jZj1]->fill(ysum(*z, j1));
      _h[kYsum2jZj2]->fill(ysum(*z, j2));
      _h[kYsum2jZdijet]->fill(ysum(*z, j1.mom() + j2.mom()));
      _h[kYsum2jJJ]->fill(ysum(j1, j2));

      fill3cuts(kDphi2jZj1, z->pt(), 150*GeV, 300*GeV, deltaPhi(*z, j1));
      _h[kDijetMass]->fill((j1.mom() + j2.mom()).mass());
      //////////////////////////////

      ////////////////////////////////
      /// Nj >= 3 in |y| < 2.4
      if (goodjets24.size() < 3) return;

      const Jet& j3 = goodjets24[2];

      fill3cuts(kDphi3jZj1, z->pt(), 150*GeV, 300*GeV, deltaPhi(*z, j1));
      fill3cuts(kDphi3jZj2, z->pt(), 150*GeV, 300*GeV, deltaPhi(*z, j2));
      fill3cuts(kDphi3jZj3, z->pt(), 150*GeV, 300*GeV, deltaPhi(*z, j3));
      fill3cuts(kDphi3jJ1j2, z->pt(), 150*GeV, 300*GeV, deltaPhi(j1, j2));
      fill3cuts(kDphi3jJ1j3, z->pt(), 150*GeV, 300*GeV, deltaPhi(j1, j3));
      fill3cuts(kDphi3jJ2j3, z->pt(), 150*GeV, 300*GeV, deltaPhi(j2, j3));

      //Although not specified in the paper, the H^{jet}_{T} varible used in
      //Fig. 24 differs from H_T and is defined as the scalar sum
      //of the pt of the three leading jets
      double ht_3jets = goodjets24[0].pt() + goodjets24[1].pt() + goodjets24[2].pt();
      if(z->pt() > 150*GeV && ht_3jets > 300*GeV){
        _h[kDphi3jZj1HighHT]->fill(deltaPhi(*z, j1));
        _h[kDphi3jZj2HighHT]->fill(deltaPhi(*z, j2));
        _h[kDphi3jZj3HighHT]->fill(deltaPhi(*z, j3));
      }
      //////////////////////////////
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

      // when running in combined mode, need to average to get lepton xsec
      if (_mode == 2) norm /= 2.;

      MSG_DEBUG("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_DEBUG("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_DEBUG("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_DEBUG("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      unsigned ih = 0;
      for(auto& h: _h){
        if(ih < nNormalized) normalize(h);
        else scale(h, norm);
        ++ih;
      }

    }

  protected:

    size_t _mode;

  private:

    ///  Histograms
    std::vector<Histo1DPtr> _h;

    /// List of histogram in the order of appearance
    /// in the paper. Histograms are numbered according
    /// to this order.
    std::vector<enum histIds> _histListInPaperOrder;

    //@name Binning of pseudo-2D/3D histograms
    std::vector<double> _ptjYjBins;
    std::vector<double> _yJyZbins;
    std::vector<double> _ptJyJyZbinsYj;
    std::vector<double>  _ptJyJyZbinsYZ;

  };

  RIVET_DECLARE_PLUGIN(CMS_2017_I1497519);
}

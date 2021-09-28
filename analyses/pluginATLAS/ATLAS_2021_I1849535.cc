// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {

  /// @name M4lLineshape analysis
  class ATLAS_2021_I1849535 : public Analysis {
    public:

      /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2021_I1849535);

      void init() {

        // Selection
        Cut el_fid_sel = (Cuts::abseta < 2.47) && (Cuts::pT > 7*GeV); 
        Cut mu_fid_sel = (Cuts::abseta < 2.7) && (Cuts::pT > 5*GeV);

        PromptFinalState photons(Cuts::abspid == PID::PHOTON);
        PromptFinalState elecs(Cuts::abspid == PID::ELECTRON);
        PromptFinalState muons(Cuts::abspid == PID::MUON && mu_fid_sel);
        elecs.acceptTauDecays(true);
        muons.acceptTauDecays(true);


        // Final state including all charged particles
        declare(ChargedFinalState(), "CFS");

        DressedLeptons dressed_elecs(photons, elecs, 0.1, el_fid_sel, false);
        declare(dressed_elecs, "elecs");

        declare(muons, "muons");

        // Book histos
        book(_h["m4l_paper"],      1,1,1);
        book(_h["m4l_4mu_paper"],  2,1,1);
        book(_h["m4l_4e_paper"],   3,1,1);
        book(_h["m4l_2e2mu_paper"],4,1,1);

        book(_h["mZ1_Z_paper"],5,1,1);  
        book(_h["mZ1_H_paper"],6,1,1);
        book(_h["mZ1_offshell_paper"],7,1,1);
        book(_h["mZ1_ZZ_paper"],8,1,1);

        book(_h["mZ2_Z_paper"],9,1,1);
        book(_h["mZ2_H_paper"],10,1,1);
        book(_h["mZ2_offshell_paper"],11,1,1);
        book(_h["mZ2_ZZ_paper"],12,1,1);

        book(_h["ptZ1_Z_paper"],13,1,1);
        book(_h["ptZ1_H_paper"],14,1,1);
        book(_h["ptZ1_offshell_paper"],15,1,1); 
        book(_h["ptZ1_ZZ_paper"],16,1,1);

        book(_h["ptZ2_Z_paper"],17,1,1);
        book(_h["ptZ2_H_paper"],18,1,1);
        book(_h["ptZ2_offshell_paper"],19,1,1);
        book(_h["ptZ2_ZZ_paper"],20,1,1);

        book(_h["costhetastar1_Z_paper"],21,1,1);
        book(_h["costhetastar1_H_paper"],22,1,1);
        book(_h["costhetastar1_offshell_paper"],23,1,1);
        book(_h["costhetastar1_ZZ_paper"],24,1,1);

        book(_h["costhetastar2_Z_paper"],25,1,1);
        book(_h["costhetastar2_H_paper"],26,1,1);
        book(_h["costhetastar2_offshell_paper"],27,1,1);
        book(_h["costhetastar2_ZZ_paper"],28,1,1);

        book(_h["dy_Z1Z2_Z_paper"]  ,29,1,1);
        book(_h["dy_Z1Z2_H_paper"]  ,30,1,1);
        book(_h["dy_Z1Z2_offshell_paper"],31,1,1);
        book(_h["dy_Z1Z2_ZZ_paper"]  ,32,1,1);

        book(_h["dphi_Z1Z2_Z_paper"]  ,33,1,1);
        book(_h["dphi_Z1Z2_H_paper"]  ,34,1,1);
        book(_h["dphi_Z1Z2_offshell_paper"],35,1,1);
        book(_h["dphi_Z1Z2_ZZ_paper"]  ,36,1,1);

        book(_h["dphi_l1l2_Z_paper"]   ,37,1,1);
        book(_h["dphi_l1l2_H_paper"]  ,38,1,1);
        book(_h["dphi_l1l2_offshell_paper"],39,1,1);
        book(_h["dphi_l1l2_ZZ_paper"]  ,40,1,1);

        book(_h["m4l_ptslice1_paper"],41,1,1);
        book(_h["m4l_ptslice2_paper"],42,1,1);
        book(_h["m4l_ptslice3_paper"],43,1,1);
        book(_h["m4l_ptslice4_paper"],44,1,1);
        book(_h["m4l_ptslice5_paper"],45,1,1);

        book(_h["m4l_yslice1_paper"],46,1,1);
        book(_h["m4l_yslice2_paper"],47,1,1);
        book(_h["m4l_yslice3_paper"],48,1,1);
        book(_h["m4l_yslice4_paper"],49,1,1);
        book(_h["m4l_yslice5_paper"],50,1,1);

      }

      /// Generic dilepton candidate
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
      bool passesTruthIsolation(Quadruplet quad, const Particles charged_tracks, Particles& truthLeptons ){
        bool pass =true;
        Particles leps;
        leps.push_back(quad._z1.first);
        leps.push_back(quad._z2.first);
        leps.push_back(quad._z1.second);
        leps.push_back(quad._z2.second);
        for (auto &lep : leps){
          double pTinCone = -lep.pT();
          for (const Particle& track : charged_tracks) {
            if (deltaR(lep.momentum(), track.momentum()) < 0.3)
              pTinCone += track.pT();
          }
          for (const Particle& tlep: truthLeptons) {
            float dR= deltaR(lep.momentum(),  tlep.momentum()); 
            if ( dR>0 && dR < 0.3)
              pTinCone -= tlep.pT();
          }
          if (pTinCone > 0.16* lep.pT()){
            pass=false;
          }
        }
        return pass;
      }

      std::vector<Quadruplet> getBestQuads(Particles& particles, bool drcut = true) {
        // H->ZZ->4l pairing
        // - Two same flavor opposite charged leptons
        // - Ambiguities in pairing are resolved by choosing the combination
        //     that results in the smaller value of |mll - mZ| for each pair successively
        std::vector<Quadruplet> quads {};

        size_t n_parts = particles.size();
        if (n_parts < 4)  return quads; 

        // STEP 1: find SFOS pairs 
        std::vector<Dilepton> SFOS;
        for (size_t i = 0; i < n_parts; ++i) {
          for (size_t j = 0; j < i; ++j) {
            if (particles[i].pid() == -particles[j].pid()) {
              // sort such that the negative lepton is listed first
              Dilepton sfos;
              if (particles[i].pid() > 0)  sfos = make_pair(particles[i], particles[j]);
              else                         sfos = make_pair(particles[j], particles[i]);

              if (sfos.mom().mass() > _ll_mass*GeV && (!drcut || deltaR(particles[i],particles[j]) > _dRll)) SFOS.push_back(sfos);
            }
          }
        }
        if (SFOS.size() < 2)  return quads;

        // now we sort the SFOS pairs
        std::sort(SFOS.begin(), SFOS.end(), [](const Dilepton& p1, const Dilepton& p2) {
            return fabs(p1.mom().mass() - Z_mass) < fabs(p2.mom().mass() - Z_mass);
            });

        //form all possible quadruplets, passing the pt cuts, the dR cuts and the mll cuts
        for (size_t k = 0; k < SFOS.size(); ++k) {
          for (size_t l = k+1; l < SFOS.size(); ++l) {
            if(drcut) {
              if (deltaR(SFOS[k].first.mom(),  SFOS[l].first.mom())  < _dRll)  continue;
              if (deltaR(SFOS[k].first.mom(),  SFOS[l].second.mom()) < _dRll)  continue;
              if (deltaR(SFOS[k].second.mom(), SFOS[l].first.mom())  < _dRll)  continue;
              if (deltaR(SFOS[k].second.mom(), SFOS[l].second.mom()) < _dRll)  continue;
            }
            if ( (SFOS[k].first.pid()   == -SFOS[l].first.pid())  && ((SFOS[k].first.mom()  + SFOS[l].first.mom()).mass()  < _ll_mass*GeV)) continue;
            if ( (SFOS[k].first.pid()   == -SFOS[l].second.pid()) && ((SFOS[k].first.mom()  + SFOS[l].second.mom()).mass() < _ll_mass*GeV)) continue;
            if ( (SFOS[k].second.pid()  == -SFOS[l].first.pid())  && ((SFOS[k].second.mom() + SFOS[l].first.mom()).mass()  < _ll_mass*GeV)) continue;
            if ( (SFOS[k].second.pid()  == -SFOS[l].second.pid()) && ((SFOS[k].second.mom() + SFOS[l].second.mom()).mass() < _ll_mass*GeV)) continue;

            //think technically this should happen before quad formation now and with all leptons not just those in quad so commenting out
            //std::vector<double> lep_pt { SFOS[k].pTl1(), SFOS[k].pTl2(), SFOS[l].pTl1(), SFOS[l].pTl2() };
            //std::sort(lep_pt.begin(), lep_pt.end(), std::greater<double>()); 
            //if (!(lep_pt[0] > _pt_lep1*GeV && lep_pt[1] > _pt_lep2*GeV && lep_pt[2] > _pt_lep3*GeV)) continue;
            quads.push_back( Quadruplet(SFOS[k], SFOS[l]) );
          }
        }
        return quads;
      }

      bool passPtLeptons(const Particles& particles) {
        size_t n_parts = particles.size();
        if (n_parts < 4)  return false;
        // cut on pT of leptons
        return ( particles[0].mom().pt() > _pt_lep1*GeV && particles[1].mom().pt() > _pt_lep2*GeV && particles[2].mom().pt() > _pt_lep3*GeV) ;
      }


      // Do the analysis
      void analyze(const Event& event) {

        const Particles charged_tracks    = apply<ChargedFinalState>(event, "CFS").particles();


        //preselection of leptons for ZZ-> llll final state
        Particles dressed_leptons;
        for (auto lep : apply<FinalState>(event, "muons").particles()) { dressed_leptons.push_back(lep); }
        for (auto lep : apply<DressedLeptons>(event, "elecs").dressedLeptons()) { dressed_leptons.push_back(lep); }



        // sort to put highest pT first
        std::sort(dressed_leptons.begin(), dressed_leptons.end(), [](const Particle& l1, const Particle& l2) {
            return l1.pt() > l2.pt();
            });

        auto foundDressedNoDrll = getBestQuads(dressed_leptons,false);

        //now doing pt cut before quad formation so also apply this here
        if (!passPtLeptons(dressed_leptons)) vetoEvent;

        auto foundDressed = getBestQuads(dressed_leptons);
        // if we don't find any quad, we can stop here 
        if (foundDressed.empty())  vetoEvent;
        if (!passesTruthIsolation(foundDressed[0], charged_tracks, dressed_leptons)) vetoEvent;

        double m4l = foundDressed[0].mom().mass()/GeV;
        double pt4l = foundDressed[0].mom().pT()/GeV;
        double y4l = foundDressed[0].mom().absrap();
        double mZ1 = foundDressed[0].getZ1().mom().mass()/GeV;
        double mZ2 = foundDressed[0].getZ2().mom().mass()/GeV;
        double ptZ1 = foundDressed[0].getZ1().mom().pT()/GeV;
        double ptZ2 = foundDressed[0].getZ2().mom().pT()/GeV;
        double dy_Z1Z2 = fabs(foundDressed[0].getZ1().mom().rapidity() - foundDressed[0].getZ2().mom().rapidity());
        double dphi_Z1Z2 = deltaPhi(foundDressed[0].getZ1().mom(),foundDressed[0].getZ2().mom());
        double dphi_l1l2 = deltaPhi(dressed_leptons[0].mom(),dressed_leptons[1].mom());

        _h["m4l_paper"]->fill(m4l);
        if (     pt4l <  10.)	  _h["m4l_ptslice1_paper"]->fill(m4l);
        else if (pt4l <  20.)	  _h["m4l_ptslice2_paper"]->fill(m4l);
        else if (pt4l < 50.)	  _h["m4l_ptslice3_paper"]->fill(m4l);
        else if (pt4l < 100.)	  _h["m4l_ptslice4_paper"]->fill(m4l);
        else if (pt4l < 600.)	  _h["m4l_ptslice5_paper"]->fill(m4l); 

        if (y4l < 0.3)	  _h["m4l_yslice1_paper"]->fill(m4l);
        else if (y4l < 0.6)	  _h["m4l_yslice2_paper"]->fill(m4l);
        else if (y4l < 0.9)	  _h["m4l_yslice3_paper"]->fill(m4l);
        else if (y4l < 1.2)	  _h["m4l_yslice4_paper"]->fill(m4l);
        else if (y4l < 2.5)	  _h["m4l_yslice5_paper"]->fill(m4l);


        Quadruplet::FlavCombi flavour = foundDressed[0].type();
        if (     flavour == Quadruplet::FlavCombi::mm) { _h["m4l_4mu_paper"]->fill(m4l); }
        else if (flavour == Quadruplet::FlavCombi::ee) { _h["m4l_4e_paper"]->fill(m4l); }
        else if (flavour == Quadruplet::FlavCombi::me || flavour == Quadruplet::FlavCombi::em) {
          _h["m4l_2e2mu_paper"]->fill(m4l);
        }

        // polarization variables
        // Get four-momentum of the first lepton pair
        const FourMomentum pcom = foundDressed.at(0).getZ1().mom();
        const Vector3 betacom = pcom.betaVec();
        const Vector3 unitboostvec = betacom.unit();
        const LorentzTransform comboost = LorentzTransform::mkFrameTransformFromBeta(betacom);
        // Get four-momentum of the negative lepton w.r.t. the first lepton pair
        const FourMomentum p1com = comboost.transform(foundDressed.at(0).getZ1().first.mom());
        float costhetastar1 = cos(p1com.p3().angle(unitboostvec));

        // Get four-momentum of the second lepton pair
        const FourMomentum pcom2 = foundDressed.at(0).getZ2().mom();
        const Vector3 betacom2 = pcom2.betaVec();
        const Vector3 unitboostvec2 = betacom2.unit();
        const LorentzTransform comboost2 = LorentzTransform::mkFrameTransformFromBeta(betacom2);
        // Get four-momentum of the negative lepton w.r.t. the second lepton pair
        const FourMomentum p2com = comboost2.transform(foundDressed.at(0).getZ2().first.mom());
        float  costhetastar2 = cos(p2com.p3().angle(unitboostvec2));

        //fill m4l binned variables
        if (60 < m4l && m4l < 100.) {
          _h["mZ1_Z_paper"]->fill(mZ1);
          _h["mZ2_Z_paper"]->fill(mZ2);
          _h["ptZ1_Z_paper"]->fill(ptZ1);
          _h["ptZ2_Z_paper"]->fill(ptZ2);
          _h["dy_Z1Z2_Z_paper"]->fill(dy_Z1Z2);
          _h["dphi_Z1Z2_Z_paper"]->fill(dphi_Z1Z2);
          _h["dphi_l1l2_Z_paper"]->fill(dphi_l1l2);
          _h["costhetastar1_Z_paper"]->fill(costhetastar1 );
          _h["costhetastar2_Z_paper"]->fill(costhetastar2 );
        }
        else if(120 < m4l && m4l < 130 ){
          _h["mZ1_H_paper"]->fill(mZ1);
          _h["mZ2_H_paper"]->fill(mZ2);
          _h["ptZ1_H_paper"]->fill(ptZ1);
          _h["ptZ2_H_paper"]->fill(ptZ2);
          _h["dy_Z1Z2_H_paper"]->fill(dy_Z1Z2);
          _h["dphi_Z1Z2_H_paper"]->fill(dphi_Z1Z2);
          _h["dphi_l1l2_H_paper"]->fill(dphi_l1l2);
          _h["costhetastar1_H_paper"]->fill(costhetastar1 );
          _h["costhetastar2_H_paper"]->fill(costhetastar2 );
        }
        else if(180 < m4l && m4l < 2000){
          _h["mZ1_ZZ_paper"]->fill(mZ1);
          _h["mZ2_ZZ_paper"]->fill(mZ2);
          _h["ptZ1_ZZ_paper"]->fill(ptZ1);
          _h["ptZ2_ZZ_paper"]->fill(ptZ2);
          _h["dy_Z1Z2_ZZ_paper"]->fill(dy_Z1Z2);
          _h["dphi_Z1Z2_ZZ_paper"]->fill(dphi_Z1Z2);
          _h["dphi_l1l2_ZZ_paper"]->fill(dphi_l1l2);
          _h["costhetastar1_ZZ_paper"]->fill(costhetastar1 );
          _h["costhetastar2_ZZ_paper"]->fill(costhetastar2 );
        }
        else{
          _h["mZ1_offshell_paper"]->fill(mZ1);
          _h["mZ2_offshell_paper"]->fill(mZ2);
          _h["ptZ1_offshell_paper"]->fill(ptZ1);
          _h["ptZ2_offshell_paper"]->fill(ptZ2);
          _h["dy_Z1Z2_offshell_paper"]->fill(dy_Z1Z2);
          _h["dphi_Z1Z2_offshell_paper"]->fill(dphi_Z1Z2);
          _h["dphi_l1l2_offshell_paper"]->fill(dphi_l1l2);
          _h["costhetastar1_offshell_paper"]->fill(costhetastar1 );
          _h["costhetastar2_offshell_paper"]->fill(costhetastar2 );
        }

      }//end analysis

      /// Finalize
      void finalize() {
        const double sf = crossSection() / femtobarn / sumOfWeights();
        for (auto hist : _h) { scale(hist.second, sf); }
      }

    private:

      map<string, Histo1DPtr> _h;
      static constexpr double Z_mass = 91.1876;
      static constexpr float _pt_lep1 = 20.;
      static constexpr float _pt_lep2 = 10.;
      static constexpr float _pt_lep3 = 0.;
      static constexpr float _ll_mass = 5.;
      static constexpr float _dRll = 0.05;

  };  // end class ATLAS_2021_I1849535

  DECLARE_RIVET_PLUGIN(ATLAS_2021_I1849535);
}  // end namespace rivet

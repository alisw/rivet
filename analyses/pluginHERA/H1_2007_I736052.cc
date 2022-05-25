// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of inclusive production of D* mesons both with and without dijet production in DIS collisions (H1)
  class H1_2007_I736052 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2007_I736052);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
        
        declare(DISLepton(), "Lepton");
        declare(DISKinematics(), "Kinematics");      
      
        const FinalState fs;
        declare(fs, "FS");
        const UnstableParticles ufs;
        declare(ufs, "UFS");
      
        double jet_radius = 1.0;
        const DISFinalState DISfs(DISFinalState::BoostFrame::BREIT);
        declare(FastJets(DISfs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::Et_scheme, jet_radius), "jets_fs");
      
        Histo1DPtr tmp;
        book(_h["411"], 4, 1, 1);
        book(_h["511"], 5, 1, 1);
        book(_h["611"], 6, 1, 1);
        book(_h["711"], 7, 1, 1);
        book(_h["811"], 8, 1, 1);
        book(_h["911"], 9, 1, 1);
        _h_binned["Q2xbj"].add( 2., 4.22, book(_h["1011"], 10, 1, 1));        
        _h_binned["Q2xbj"].add( 4.22, 10., book(_h["1111"], 11, 1, 1));
        _h_binned["Q2xbj"].add( 10., 17.8, book(_h["1211"], 12, 1, 1));
        _h_binned["Q2xbj"].add( 17.8, 31.6, book(_h["1311"], 13, 1, 1));
        _h_binned["Q2xbj"].add( 31.6, 100., book(_h["1411"], 14, 1, 1));
        book(_h["1511"], 15, 1, 1);
        book(_h["1611"], 16, 1, 1);
        book(_h["1711"], 17, 1, 1);
        book(_h["1811"], 18, 1, 1);
        book(_h["1911"], 19, 1, 1);
        book(_h["2011"], 20, 1, 1);
        book(_h["2111"], 21, 1, 1);
        book(_h["2211"], 22, 1, 1);

        _h_binned["Q2phi"].add( 2., 10., book(tmp, 23, 1, 1));
        _h_binned["Q2phi"].add( 10., 100., book(tmp, 24, 1, 1));
        book(_h["2511"], 25, 1, 1);
        book(_h["2611"], 26, 1, 1);
        book(_h["2711"], 27, 1, 1);
        book(_h["2811"], 28, 1, 1);
        _h_binned["Q2xgam"].add( 2., 5., book(_h["2911"], 29, 1, 1));
        _h_binned["Q2xgam"].add( 5., 10., book(_h["2912"], 29, 1, 2));
        _h_binned["Q2xgam"].add( 10., 100., book(_h["2913"], 29, 1, 3));
        book(_h["3011"], 30, 1, 1);
        _h_binned["Q2xglue"].add( 2., 5., book(_h["3111"], 31, 1, 1));
        _h_binned["Q2xglue"].add( 5., 10., book(_h["3112"], 31, 1, 2));
        _h_binned["Q2xglue"].add( 10., 100., book(_h["3113"], 31, 1, 3));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        
        const FinalState& fs = apply<FinalState>(event, "FS");
        const size_t numParticles = fs.particles().size();
        
        Jets jets_fs = apply<JetAlg>(event, "jets_fs").jetsByPt(); // Jets with cut on eta
        double jet_radius = 1.0;
        
        const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
        
        if (numParticles < 2){
            MSG_DEBUG("Failed leptonic event cut");
            vetoEvent;
        }
        
        Particles Dstar;
        for(const Particle& p : filter_select(ufs.particles(), Cuts::pT > 1.5*GeV and Cuts::pT < 15*GeV and Cuts::abseta < 1.5 and Cuts::abspid==413)) {
            Dstar.push_back(p);
        }
        if(Dstar.size() == 0){ // Cut on Dstar
            MSG_DEBUG("Failed Dstar cut");
            vetoEvent;
        }
            
        const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
        const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");
        
       
        double Q2 = dk.Q2();
        double y  = dk.y();

        if(y < 0.05 or y > 0.7 or Q2 < 2 or Q2 > 100){ // Cut on event kinematics
            MSG_DEBUG("Failed kinematics cut");
            vetoEvent; 
        } 
        
        // Extract the particles other than the lepton
        Particles particles;
        particles.reserve(fs.particles().size());
        ConstGenParticlePtr dislepGP = dl.out().genParticle();
        for (const Particle& p : fs.particles()) {
           ConstGenParticlePtr loopGP = p.genParticle();
           if (loopGP == dislepGP) continue;
           particles.push_back(p);
        }

        
        
        const LorentzTransform hcmboost = dk.boostHCM(); // Hadron cm system 
        const LorentzTransform breitboost = dk.boostBreit(); //Breit system
        const LorentzTransform labboost = breitboost.inverse(); //Labsystem from Breit
        
        double xbj  = dk.x();
        double W2 = dk.W2();
        double pT_cm(0); 
        
        _h["411"] -> fill(Q2);
        _h["511"] -> fill(xbj);
        _h["611"] -> fill(std::sqrt(W2)); 
        
        _h_binned["Q2xbj"].fill(Q2, xbj);      
        
        for(const Particle& p : Dstar){
            _h["711"] -> fill(p.pT());
            _h["811"] -> fill(p.eta());
            _h["911"] -> fill((p.E() - p.pz())/(2*y*dk.beamLepton().E()));
            
            const FourMomentum hcmMom = hcmboost.transform(p.momentum()); 
            if(pT_cm < hcmMom.pT()) pT_cm = hcmMom.pT();
            
            if(hcmMom.pT() > 2){              //Cut for 1711 and 1811 histograms
                _h["1711"] -> fill(p.pT());
                _h["1811"] -> fill(p.eta());
            }            
        }

        if(pT_cm > 2){
            _h["1511"] -> fill(Q2);
            _h["1611"] -> fill(xbj);
        }     
        
/*
        FourMomentum gammaZ;
        const FourMomentum leptonOUT = dl.out();
        const FourMomentum leptonIN = dl.in();
        gammaZ = leptonIN - leptonOUT ;
        cout << " LAB: gammaZ " <<  gammaZ <<  endl; 
        cout << " HCM: gammaZ " <<   hcmboost.transform(gammaZ) <<  endl; 
        cout << " LAB: p-beam " <<  dk.beamHadron().momentum() <<  endl; 
        cout << " HCM: p-beam " <<   hcmboost.transform(dk.beamHadron().momentum()) <<  endl; 
        cout << " Breit: gammaZ " <<   breitboost.transform(gammaZ) <<  endl; 
        cout << " Breit: p-beam " <<   breitboost.transform(dk.beamHadron().momentum()) <<  endl; 
        for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())) {
          const PdgId pid = p->pdg_id();
          if (abs(pid) == 23) {
          cout<< " HEPMC: photon/Z " << p->momentum() << endl;
          }
        }        
*/
        bool two_jets_with_cut(false); 
        
        Jets jet_cut; 
        
        for(const Jet& j : jets_fs){
            double etalab = labboost.transform(j.momentum()).eta() ;
            if(j.momentum().Et() > 3 and etalab > -1 and etalab < 2.5){  // Cut on jet energy in Breit sys.
                jet_cut.push_back(j);                            }     
        }  

        if(jet_cut.size() >= 2 && jet_cut[0].momentum().Et() > 4 ) two_jets_with_cut = true;
        
        double delta_phi;
        FourMomentum jet1;
        FourMomentum jet2;
        
        Jets jet_DJ, jet_OJ;
        bool found_DJ(false), found_OJ(false);
        //cout << " check D and other jets " << endl;        
        if(two_jets_with_cut){
            jet1 = jet_cut[0].momentum(); // momentum of jet #1 in Breit sys. 
            jet2 = jet_cut[1].momentum(); // momentum of jet #2 in Breit sys. 

            delta_phi = deltaPhi(jet1,jet2)/degree ;
            _h["1911"] -> fill(Q2);     
            _h["2011"] -> fill(xbj); 
            _h["2111"] -> fill(jet_cut[0].momentum().Et()/GeV); 
            _h["2211"] -> fill(FourMomentum(jet1+jet2).mass()/GeV);  // Jets invariant mass in Breit sys.
            //cout << " dphi " << delta_phi <<" " << deltaPhi(jet_cut[0].momentum(),jet_cut[1].momentum())/degree <<  endl;
            _h_binned["Q2phi"].fill(Q2, delta_phi);
            for (const Jet& jet : jet_cut) {
                for(const Particle & p : Dstar) {
                    if(deltaR(breitboost.transform(p.momentum()), jet.momentum()) < jet_radius  ) {
                        jet_DJ.push_back(jet);
                        found_DJ = true;
                    }
                    else {
                    jet_OJ.push_back(jet);
                    found_OJ = true;
                    }
                }                

               if( found_DJ &&  found_OJ ) break ;  
            }
            //cout << " DJ jet size " << jet_DJ.size() << "found Dj " << found_DJ <<  " OJ jet size " << jet_DJ.size()  << " found OJ " << found_OJ << endl;
            if(jet_DJ.size()>0 && jet_OJ.size()>0 )
            {            
                double eta_DJ_breit = jet_DJ[0].momentum().eta();
                double eta_OJ_breit = jet_OJ[0].momentum().eta(); 
                
                _h["2511"] -> fill(eta_DJ_breit);
                _h["2611"] -> fill(eta_OJ_breit);
                _h["2711"] -> fill(abs(eta_DJ_breit - eta_OJ_breit));    
                
                double x_gamma, x_gluon;
                double E_star_p_z_star(0); // E() - pz() sum for for x_gamma in \gamma p cm
                double E_star_p_z_had(0);

                // for jets in had cms: boost first back to lab and then to hcm
               Jet jet_DJ_hcm, jet_OJ_hcm;
               jet_DJ_hcm = hcmboost.transform(labboost(jet_DJ[0].momentum()));
               jet_OJ_hcm = hcmboost.transform(labboost(jet_OJ[0].momentum()));

// Need to change sign: by default hcmboost has gamma* in +z dir, and p in -z dir, but here we need: gamma* -in -z and proton in +z.
// so - -> +: watch out for E-pz -> E+pz  and exp(eta_OJ_hcm) -> exp(-eta_OJ_hcm)
// this was already noted in  H1_2007_I746380.cc


   
               double Et_DJ_hcm = jet_DJ_hcm.Et();
               double eta_DJ_hcm = jet_DJ_hcm.eta();
               double Et_OJ_hcm = jet_OJ_hcm.Et();
               double eta_OJ_hcm = jet_OJ_hcm.eta();

                //observed fraction of the photon momentum carried by the parton involved in the hard subprocess 
                
                E_star_p_z_star = jet_DJ_hcm.momentum().E() + jet_DJ_hcm.momentum().pz() + jet_OJ_hcm.momentum().E() + jet_OJ_hcm.momentum().pz();

                for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
                  const Particle& p = particles[ip1];
                  E_star_p_z_had += (hcmboost.transform(p.momentum())).E() + (hcmboost.transform(p.momentum())).pz();
                }
                x_gamma = E_star_p_z_star/E_star_p_z_had;
                x_gluon = (Et_OJ_hcm*exp(-eta_OJ_hcm) + Et_DJ_hcm*exp(-eta_DJ_hcm))/(2.*hcmboost.transform(dk.beamHadron().momentum()).E());
                //observed fraction of the proton momentum carried by the gluon 
                
                //cout << " xgamm " << x_gamma << " x_glu " << x_gluon << endl;
                _h["2811"] -> fill(x_gamma); 
                _h_binned["Q2xgam"].fill(Q2, x_gamma);  
                
                _h["3011"] -> fill(log10(x_gluon));
                _h_binned["Q2xglue"].fill(Q2, log10(x_gluon));                  
            }
        }        
    }


    /// Normalise histograms etc., after the run
    void finalize() {
        double norm = crossSection()/nanobarn/sumW(); 
        //double norm = crossSection()/nanobarn/sumOfWeights(); 
        //cout << " SumOfWeights " << sumW() << " "<< sumOfWeights() << endl;
       
        scale(_h["411"], norm);
        scale(_h["511"], norm);
        scale(_h["611"], norm);
        scale(_h["711"], norm);
        scale(_h["811"], norm);
        scale(_h["911"], norm);
        _h_binned["Q2xbj"].scale(norm, this);
        scale(_h["1511"], norm);
        scale(_h["1611"], norm);
        scale(_h["1711"], norm);
        scale(_h["1811"], norm);
        scale(_h["1911"], norm);
        scale(_h["2011"], norm);
        scale(_h["2111"], norm);
        scale(_h["2211"], norm);
        _h_binned["Q2phi"].scale(norm*180./M_PI, this);
        scale(_h["2511"], norm);
        scale(_h["2611"], norm);
        scale(_h["2711"], norm);
        scale(_h["2811"], norm);
        _h_binned["Q2xgam"].scale(norm, this);
        
        scale(_h["3011"], norm);
        // new scaling needed, since x bins are in log10(x)
        vector<YODA::HistoBin1D>& bins = _h["3011"] -> bins(); 
        for (auto & b : bins) {          
           double scale_new = b.xWidth()/pow(10,b.xWidth()) ;
           // jacobian d log10/dx = dlog10 / dlogx dlogx/dx = 2.3026 /x 
           double factor = pow(10,b.xMid())/2.3026;
           b.scaleW(scale_new/factor) ;
        }
        
        _h_binned["Q2xglue"].scale(norm, this);
        // new scaling needed, since x bins are in log10(x)
        for (Histo1DPtr histo : _h_binned["Q2xglue"].histos()) {
           vector<YODA::HistoBin1D>& bins = histo -> bins(); 
           for (auto & b : bins) {             
              double scale_new = b.xWidth()/pow(10,b.xWidth()) ;
              // jacobian d log10/dx = dlog10 / dlogx dlogx/dx = 2.3026 /x 
              double factor = pow(10,b.xMid())/2.3026;
              b.scaleW(scale_new/factor) ;
           }

        }

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, BinnedHistogram> _h_binned;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(H1_2007_I736052);

}

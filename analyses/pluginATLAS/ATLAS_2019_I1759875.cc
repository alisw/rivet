#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  /// @brief lepton differential ttbar analysis at 13 TeV
  class ATLAS_2019_I1759875 : public Analysis {
  public:
    
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1759875);
    
    void init() {

      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;

      // All final state particles
      const FinalState fs;

      // Get photons to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      PromptFinalState prompt_el(Cuts::abspid == PID::ELECTRON, true);
      DressedLeptons elecs(photons, prompt_el, 0.1, (Cuts::abseta < 2.5) && (Cuts::pT > 20*GeV));
      DressedLeptons veto_elecs(photons, prompt_el, 0.1, eta_full, false);
      declare(elecs, "elecs");

      // Projection to find the muons
      PromptFinalState prompt_mu(Cuts::abspid == PID::MUON, true);
      DressedLeptons muons(photons, prompt_mu, 0.1, (Cuts::abseta < 2.5) && (Cuts::pT > 20*GeV));
      DressedLeptons veto_muons(photons, prompt_mu, 0.1, eta_full, false);
      declare(muons, "muons");

      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(veto_elecs);
      vfs.addVetoOnThisFinalState(veto_muons);

      // Book histograms
      bookHistos("lep_pt",       1);
      bookHistos("lep_eta",      3);
      bookHistos("dilep_pt",     5);
      bookHistos("dilep_mass",   7);
      bookHistos("dilep_rap",    9);
      bookHistos("dilep_dphi",  11);
      bookHistos("dilep_sumpt", 13);
      bookHistos("dilep_sumE",  15);

      // unrolled 2D distributions - 2nd-dim bin edges must be specified
      std::vector<double> massbins={0.,80.,120.,200.,500.};
      
      bookHisto2D("lep_eta_mass",17,massbins);
      bookHisto2D("dilep_rap_mass",19,massbins);
      bookHisto2D("dilep_dphi_mass",21,massbins);
    }

    void analyze(const Event& event) {
      vector<DressedLepton> elecs = apply<DressedLeptons>(event, "elecs").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "muons").dressedLeptons();

      if (elecs.empty() || muons.empty())  vetoEvent;
      if (elecs[0].charge() == muons[0].charge())  vetoEvent;  
 
      FourMomentum el = elecs[0].momentum();
      FourMomentum mu = muons[0].momentum();
      FourMomentum ll = elecs[0].momentum() + muons[0].momentum();
                  
      // Fill histograms
      // include explicit overflow protection as last bins are inclusive
      fillHistos("lep_pt",      min(el.pT()/GeV,299.));
      fillHistos("lep_pt",      min(mu.pT()/GeV,299.));
      fillHistos("lep_eta",     el.abseta());
      fillHistos("lep_eta",     mu.abseta());
      fillHistos("dilep_pt",    min(ll.pT()/GeV,299.));
      fillHistos("dilep_mass",  min(ll.mass()/GeV,499.));
      fillHistos("dilep_rap",   ll.absrap());
      fillHistos("dilep_dphi",  deltaPhi(el, mu));
      fillHistos("dilep_sumpt", min((el.pT()+mu.pT())/GeV,399.));
      fillHistos("dilep_sumE",  min((el.E()+mu.E())/GeV,699.));

      // find mass bin variable, with overflow protection
      float massv=ll.mass()/GeV;
      if (massv>499.) massv=499.;
      // Fill unrolled 2D histograms vs mass
      fillHisto2D("lep_eta_mass",el.abseta(),massv);
      fillHisto2D("lep_eta_mass",mu.abseta(),massv);
      fillHisto2D("dilep_rap_mass",ll.absrap(),massv);
      fillHisto2D("dilep_dphi_mass",deltaPhi(el,mu),massv);
    }

    void finalize() {
      // Normalize to cross-section
      const double sf = crossSection()/femtobarn/sumOfWeights();

      // finalisation of 1D histograms
      for (auto& hist : _h) {
        const double norm = 1.0 / hist.second->integral();
        // histogram normalisation
        if (hist.first.find("norm") != string::npos)  scale(hist.second, norm);
        else  scale(hist.second, sf);
      }

      // finalisation of 2D histograms
      for (auto& hist : _h_multi) {
        if (hist.first.find("_norm") != std::string::npos) {
          // scaling for normalised distribution according integral of whole set
          double norm2D = integral2D(hist.second);
          hist.second.scale(1./norm2D, this);
        }
        else {
          // scaling for non-normalised distribution
          hist.second.scale(sf, this);
        }
      }
    }


  private:

    /// @name Histogram helper functions
    //@{
    void bookHistos(const std::string name, unsigned int index) {
      book(_h[name], index, 1, 1);
      book(_h["norm_" + name],index + 1, 1, 1);
    }

    void fillHistos(const std::string name, double value) {
      _h[name]->fill(value);
      _h["norm_" + name]->fill(value);
    }

    void bookHisto2D(const std::string name, unsigned int index, std::vector<double> massbins) {
      unsigned int nmbins = massbins.size()-1;
      for (unsigned int i=0; i < nmbins; ++i) {
	{ Histo1DPtr tmp; _h_multi[name].add(massbins[i], massbins[i+1], book(tmp, index,1,1+i)); }
	const std::string namen = name+"_norm";
	{ Histo1DPtr tmp; _h_multi[namen].add(massbins[i], massbins[i+1], book(tmp, index+1,1,1+i)); }
      }
    }



    void fillHisto2D(const std::string name,
		     double val, double massval) {
      _h_multi[name].fill(massval, val);
      _h_multi[name+"_norm"].fill(massval, val);
    }


    double integral2D(BinnedHistogram& h_multi) {
        double total_integral = 0;
        for  (Histo1DPtr& h : h_multi.histos()) {
          total_integral += h->integral(false);
        }
        return total_integral;
      }


    // pointers to 1D and 2D histograms
    map<string, Histo1DPtr> _h;
    map<string, BinnedHistogram> _h_multi;
    //@}
    // acceptance counter

  };

  // Declare the class as a hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1759875);
}

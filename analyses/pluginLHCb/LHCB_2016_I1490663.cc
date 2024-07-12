// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// LHCb prompt charm hadron pT and rapidity spectra
  class LHCB_2016_I1490663 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2016_I1490663);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      declare(UnstableParticles(), "UFS");

      /// Book histograms
      /// @todo Make this interface nicer!
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add(2.0, 2.5, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add(2.5, 3.0, book(tmp, 1, 1, 2));}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add(3.0, 3.5, book(tmp, 1, 1, 3));}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add(3.5, 4.0, book(tmp, 1, 1, 4));}
      {Histo1DPtr tmp; _h_pdg411_Dplus_pT_y.add(4.0, 4.5, book(tmp, 1, 1, 5));}

      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add(2.0, 2.5, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add(2.5, 3.0, book(tmp, 2, 1, 2));}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add(3.0, 3.5, book(tmp, 2, 1, 3));}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add(3.5, 4.0, book(tmp, 2, 1, 4));}
      {Histo1DPtr tmp; _h_pdg421_Dzero_pT_y.add(4.0, 4.5, book(tmp, 2, 1, 5));}

      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add(2.0, 2.5, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add(2.5, 3.0, book(tmp, 3, 1, 2));}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add(3.0, 3.5, book(tmp, 3, 1, 3));}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add(3.5, 4.0, book(tmp, 3, 1, 4));}
      {Histo1DPtr tmp; _h_pdg431_Dsplus_pT_y.add(4.0, 4.5, book(tmp, 3, 1, 5));}

      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add(2.0, 2.5, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add(2.5, 3.0, book(tmp, 4, 1, 2));}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add(3.0, 3.5, book(tmp, 4, 1, 3));}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add(3.5, 4.0, book(tmp, 4, 1, 4));}
      {Histo1DPtr tmp; _h_pdg413_Dstarplus_pT_y.add(4.0, 4.5, book(tmp, 4, 1, 5));}

      for (int i = 0; i< 5; ++i) {
      	{Histo1DPtr tmp; _hbr_Dzero.add(2.0+i*0.5, 2.5+i*0.5, book(tmp, "TMP/Dzero_b"+to_str(i+1), refData(9, 1, 2)));}
      	{Histo1DPtr tmp; _hbr_Dplus.add(2.0+i*0.5, 2.5+i*0.5, book(tmp, "TMP/Dplus_b"+to_str(i+1), refData(9, 1, 2)));}
      	{Histo1DPtr tmp; _hbr_Ds.add(2.0+i*0.5, 2.5+i*0.5, book(tmp, "TMP/Ds_b"+to_str(i+1), refData(9, 1, 2)));}
      	{Histo1DPtr tmp; _hbr_Dstar.add(2.0+i*0.5, 2.5+i*0.5, book(tmp, "TMP/Dstar_b"+to_str(i+1), refData(9, 1, 2)));}
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Use PrimaryHadrons to avoid double counting and automatically remove the contributions from unstable?
      const UnstableParticles &ufs = apply<UnstableParticles> (event, "UFS");
      for (const Particle& p : ufs.particles() ) {

        // We're only interested in charm hadrons
        //if (!p.isHadron() || !p.hasCharm()) continue;

        PdgId apid = p.abspid();

        // do not use Cuts::abspid to avoid supplemental iteration on particles?
        if ((apid != 411) && (apid != 421) && (apid != 431) && (apid != 413)) continue;

        // Experimental selection removes non-prompt charm hadrons: we ignore those from b decays
        if (p.fromBottom()) continue;

        // Kinematic acceptance
        const double y = p.absrap(); ///< Double analysis efficiency with a "two-sided LHCb"
        const double pT = p.pT()/GeV;

        // Fiducial acceptance of the measurements
        if ((pT > 10.0) || (y < 2.0) || (y > 4.5)) continue;

        Particles daus;

        switch (apid) {
        case 411:
          _h_pdg411_Dplus_pT_y.fill(y, pT);
          // veto on decay channel [D+ -> K- pi+ pi+]cc
          if (p.children().size() != 3) break;
          if ( ((p.children(Cuts::pid == -321).size() == 1) && (p.children(Cuts::pid == 211).size() == 2)) ||
          		 ((p.children(Cuts::pid == 321).size() == 1) && (p.children(Cuts::pid == -211).size() == 2)) )
          	_hbr_Dplus.fill(y, pT); // MSG_INFO("Found [ D+ -> K- pi+ pi+ ]cc..."); };
          break;
        case 421:
          _h_pdg421_Dzero_pT_y.fill(y, pT);
          // veto on decay channel [D0 -> K- pi+]cc
          if (p.children().size() != 2) break;
          if ( ((p.children(Cuts::pid == -321).size() == 1) && (p.children(Cuts::pid == 211).size() == 1)) ||
          		 ((p.children(Cuts::pid == 321).size() == 1) && (p.children(Cuts::pid == -211).size() == 1)) )
          	_hbr_Dzero.fill(y, pT); // MSG_INFO("Found [ D0 -> K- pi+ ]cc..."); };
          break;
        case 431:
          _h_pdg431_Dsplus_pT_y.fill(y, pT);
          //veto on decay channel [Ds+ -> [K+ K-]phi0 pi+]cc
          if (p.children().size() != 2) break;
          daus = p.children(Cuts::pid == 333);
          if ( (daus.size() == 1) && (p.children(Cuts::abspid == 211).size() == 1) &&
          		 (daus.front().children(Cuts::abspid ==321).size() == 2) )
          	_hbr_Ds.fill(y, pT); // MSG_INFO("Found [ Ds+ -> phi0(-> K+ K-) pi+ ]cc..."); };
          break;
        case 413:
          _h_pdg413_Dstarplus_pT_y.fill(y, pT);
          // veto on decay channel [D*+ -> [K- pi+]D0 pi+]cc
          if (p.children().size() != 2) break;
          daus = p.children(Cuts::pid == 421);
          if ( (daus.size() == 1) && (p.children(Cuts::abspid == 211).size() == 1) &&
          		( daus.front().children().size() == 2 ) &&
          		( ( (daus.front().children(Cuts::pid == -321).size() == 1 ) && (daus.front().children(Cuts::pid == 211).size() == 1 )	) ||
          		  ( (daus.front().children(Cuts::pid == 321).size() == 1 ) && (daus.front().children(Cuts::pid == -211).size() == 1 ) ) ) )
          	_hbr_Dstar.fill(y, pT); // MSG_INFO("Found [ D*+ -> D0 (-> K- pi+)cc pi+ ]cc..."); };
          break;
        default:
        	break;
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// Factor of 0.5 to correct for the abs(rapidity) used above
      const double scale_factor = 0.5 * crossSection()/microbarn / sumOfWeights();

      /// Avoid the implicit division by the bin width in the BinnedHistogram::scale method.
      /// @todo Another thing to make nicer / more flexible in BinnedHisto
      for (Histo1DPtr h : _h_pdg411_Dplus_pT_y.histos()) h->scaleW(scale_factor);
      for (Histo1DPtr h : _h_pdg421_Dzero_pT_y.histos()) h->scaleW(scale_factor);
      for (Histo1DPtr h : _h_pdg431_Dsplus_pT_y.histos()) h->scaleW(scale_factor);
      for (Histo1DPtr h : _h_pdg413_Dstarplus_pT_y.histos()) h->scaleW(scale_factor);

      // Do ratios
      for (int i = 0; i < 5; ++i) {
      	book(hr_DplusDzero[i], 9, 1, i+1, true);
      	book(hr_DsDzero[i], 10, 1, i+1, true);
      	book(hr_DstarDzero[i], 11, 1, i+1, true);
      	book(hr_DsDplus[i], 12, 1, i+1, true);
      	book(hr_DstarDplus[i], 13, 1, i+1, true);
      	book(hr_DsDstar[i], 14, 1, i+1, true);
      	ratioScatterBins(_hbr_Dplus.histos()[i], _hbr_Dzero.histos()[i], hr_DplusDzero[i]);
      	ratioScatterBins(_hbr_Ds.histos()[i], _hbr_Dzero.histos()[i], hr_DsDzero[i]);
      	ratioScatterBins(_hbr_Dstar.histos()[i], _hbr_Dzero.histos()[i], hr_DstarDzero[i]);
      	ratioScatterBins(_hbr_Ds.histos()[i], _hbr_Dplus.histos()[i], hr_DsDplus[i]);
      	ratioScatterBins(_hbr_Dstar.histos()[i], _hbr_Dplus.histos()[i], hr_DstarDplus[i]);
      	ratioScatterBins(_hbr_Ds.histos()[i], _hbr_Dstar.histos()[i], hr_DsDstar[i]);
      	// scale 100x as measurement is in %
      	hr_DplusDzero[i]->scaleY(100.);
      	hr_DsDzero[i]->scaleY(100.);
      	hr_DstarDzero[i]->scaleY(100.);
      	hr_DsDplus[i]->scaleY(100.);
      	hr_DstarDplus[i]->scaleY(100.);
      	hr_DsDstar[i]->scaleY(100.);
      }

    }

    /// @}


  private:

    void ratioScatterBins(Histo1DPtr& hn, Histo1DPtr& hd, Scatter2DPtr &s) {
    	vector<double> sedges;
    	// extract bin edges from Scatter2D
    	for (auto p=s->points().begin(); p != s->points().end(); ++p) {
    		sedges.push_back((*p).xMin());
    		// MSG_INFO("Scatter2D bin: " << (*p).xMin() << " - " << (*p).xMax());
    	};
    	sedges.push_back(s->points().back().xMax());
    	// make deep-copies as rebinning changes bins each time - any smarter alternative ?!
    	Histo1D *hnc, *hdc;
    	hnc = new YODA::Histo1D(hn->bins(), hn->totalDbn(), hn->underflow(), hn->overflow());
    	hdc = new YODA::Histo1D(hd->bins(), hd->totalDbn(), hd->underflow(), hd->overflow());
    	hnc->rebinTo(sedges);
    	hdc->rebinTo(sedges);
    	divide(*hnc, *hdc, s);
    	delete hnc; delete hdc;
    }


    /// @name Histograms
    /// @{
    BinnedHistogram _h_pdg411_Dplus_pT_y, _hbr_Dplus;
    BinnedHistogram _h_pdg421_Dzero_pT_y, _hbr_Dzero;
    BinnedHistogram _h_pdg431_Dsplus_pT_y, _hbr_Ds;
    BinnedHistogram _h_pdg413_Dstarplus_pT_y, _hbr_Dstar;
    Scatter2DPtr hr_DplusDzero[5];
    Scatter2DPtr hr_DsDzero[5];
    Scatter2DPtr hr_DstarDzero[5];
    Scatter2DPtr hr_DsDplus[5];
    Scatter2DPtr hr_DstarDplus[5];
    Scatter2DPtr hr_DsDstar[5];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(LHCB_2016_I1490663);

}

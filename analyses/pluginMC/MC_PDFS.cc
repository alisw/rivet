// -*- C++ -*-
#include "Rivet/Analysis.hh"
// #include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  /// Generic analysis looking at various distributions of final state particles
  class MC_PDFS : public Analysis {
  public:

    /// Constructor
    MC_PDFS()
      : Analysis("MC_PDFS")
    {    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      // declare(ChargedFinalState((Cuts::etaIn(-5.0, 5.0) && Cuts::pT >=  500*MeV)), "CFS");

      // Histograms
      book(_histPdfX ,"PdfX", logspace(50, 0.000001, 1.0));
      book(_histPdfXmin ,"PdfXmin", logspace(50, 0.000001, 1.0));
      book(_histPdfXmax ,"PdfXmax", logspace(50, 0.000001, 1.0));
      book(_histPdfQ ,"PdfQ", 50, 0.0, 30.0);
      book(_histPdfXQ,"PdfXQ", logspace(50, 0.000001, 1.0), linspace(50, 0.0, 30.0));
      //book( _histPdfTrackptVsX ,"PdfTrackptVsX", logspace(50, 0.000001, 1.0));
      //book( _histPdfTrackptVsQ ,"PdfTrackptVsQ", 50, 0.0, 30.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // This analysis needs a valid HepMC PDF info object to do anything
      if (event.genEvent()->pdf_info() == 0) vetoEvent;
      PdfInfo pdfi = *(event.genEvent()->pdf_info());

#ifdef RIVET_ENABLE_HEPMC_3
      MSG_DEBUG("PDF Q = " << pdfi.scale<< " for (id, x) = "
                << "(" << pdfi.pdf_id[0] << ", " << pdfi.x[0] << ") "
                << "(" << pdfi.pdf_id[1] << ", " << pdfi.x[1] << ")");
      _histPdfX->fill(pdfi.x[0], weight);
      _histPdfX->fill(pdfi.x[1], weight);
      _histPdfXmin->fill(std::min(pdfi.x[0], pdfi.x[1]), weight);
      _histPdfXmax->fill(std::max(pdfi.x[0], pdfi.x[1]), weight);
      _histPdfQ->fill(pdfi.scale, weight); // always in GeV?
      _histPdfXQ->fill(pdfi.x[0], pdfi.scale, weight); // always in GeV?
      _histPdfXQ->fill(pdfi.x[1], pdfi.scale, weight); // always in GeV?
      
#else
      MSG_DEBUG("PDF Q = " << pdfi.scalePDF() << " for (id, x) = "
                << "(" << pdfi.id1() << ", " << pdfi.x1() << ") "
                << "(" << pdfi.id2() << ", " << pdfi.x2() << ")");
      _histPdfX->fill(pdfi.x1(), weight);
      _histPdfX->fill(pdfi.x2(), weight);
      _histPdfXmin->fill(std::min(pdfi.x1(), pdfi.x2()), weight);
      _histPdfXmax->fill(std::max(pdfi.x1(), pdfi.x2()), weight);
      _histPdfQ->fill(pdfi.scalePDF(), weight); // always in GeV?
      _histPdfXQ->fill(pdfi.x1(), pdfi.scalePDF(), weight); // always in GeV?
      _histPdfXQ->fill(pdfi.x2(), pdfi.scalePDF(), weight); // always in GeV?
#endif
      // const FinalState& cfs = apply<FinalState>(event, "CFS");
      // for (const Particle& p : cfs.particles()) {
      //   if (fabs(eta) < 2.5 && p.pT() > 10*GeV) {
      //     _histPdfTrackptVsX->fill(pdfi.x1(), p.pT()/GeV, weight);
      //     _histPdfTrackptVsX->fill(pdfi.x2(), p.pT()/GeV, weight);
      //     _histPdfTrackptVsQ->fill(pdfi.scalePDF(), p.pT()/GeV, weight);
      //   }
      // }

    }



    /// Finalize
    void finalize() {
      scale(_histPdfX, 1/sumOfWeights());
      scale(_histPdfXmin, 1/sumOfWeights());
      scale(_histPdfXmax, 1/sumOfWeights());
      scale(_histPdfQ, 1/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histPdfX, _histPdfXmin, _histPdfXmax, _histPdfQ;
    Histo2DPtr _histPdfXQ;
    // Profile1DPtr   _histPdfTrackptVsX, _histPdfTrackptVsQ;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_PDFS);

}

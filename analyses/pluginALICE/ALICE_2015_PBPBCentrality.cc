/**
 * @file   ALICE_2015_PBPBCentrality.cc
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Aug 22 15:56:23 2018
 * 
 * @brief  Dummy analysis for centrality calibration in Pb-Pb at 5.02TeV
 */
#include <Rivet/Analysis.hh>
#include <Rivet/Projections/AliceCommon.hh>

namespace Rivet
{
  /** 
   * Dummy analysis for centrality calibration in Pb-Pb at 5.02TeV
   */
  class ALICE_2015_PBPBCentrality : public Analysis
  {
  public:
    /** 
     * Constructor 
     */
    ALICE_2015_PBPBCentrality()
      : Analysis("ALICE_2015_PBPBCentrality")
    {
    }
    /** 
     * Initialize this analysis. 
     */
    void init()
    {
      ALICE::V0AndTrigger v0and;
      declare<ALICE::V0AndTrigger>(v0and,"V0-AND");

      ALICE::V0MMultiplicity v0m;
      declare<ALICE::V0MMultiplicity>(v0m,"V0M");

      _v0m = bookHisto1D("V0M","Forward multiplicity","V0M","Events");
      _imp = bookHisto1D("V0M_IMP",100,0,20,
			 "Impact parameter","b (fm)","Events");
    }
    /** 
     * Analyse a single event.
     *
     * @param event The event 
     */
    void analyze(const Event& event)
    {
      // Get and fill in the impact parameter value if the information
      // is valid.
      const HepMC::GenEvent* ge = event.genEvent();
      const HepMC::HeavyIon* hi = ge->heavy_ion();
      if (hi && hi->is_valid())
	_imp->fill(hi->impact_parameter(), event.weight());
	  

      // Check if we have any hit in either V0-A or -C.  If not, the
      // event is not selected and we get out.
      if (!apply<ALICE::V0AndTrigger>(event,"V0-AND")()) return;

      // Fill in the V0 multiplicity for this event 
      _v0m->fill(apply<ALICE::V0MMultiplicity>(event,"V0M")(), event.weight());
    }
    /** 
     * Finalize this analysis
     */
    void finalize()
    {
      _v0m->normalize();
      _imp->normalize();
    }

    /** The distribution of V0M multiplicity */
    Histo1DPtr _v0m;
    /** The distribution of impact parameters */
    Histo1DPtr _imp;
  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2015_PBPBCentrality);
}
//
// EOF
//

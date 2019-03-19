/**
 * @file   AliceCommon.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Apr 27 17:07:35 2017
 * 
 * @brief  Constants, etc. used in the ALICE Rivet code
 * @copyright GNU Lesser Public License.
 *
 * @ingroup alice_rivet
 */
#ifndef TOOLS_ALICECOMMON_HH
#define TOOLS_ALICECOMMON_HH
#include <Rivet/Tools/Cuts.hh>
#include <Rivet/Particle.hh>

namespace Rivet
{
  /** 
   * @defgroup alice_rivet ALICE specific code in Rivet 
   *
   * This include projections to emulate trigger conditions,
   * centrality, and selection of primary particles. 
   */
  /** 
   * Namespace for ALICE specific core code 
   *
   * @ingroup alice_rivet
   */
  namespace ALICE
  {
    /**
     * The acceptance cut for the V0A
     *
     * @ingroup alice_rivet
     */
    const Cut V0Aacceptance = (Cuts::etaIn(+2.8,+5.1)&&(Cuts::abscharge3 > 0));
    /**
     * The acceptance cut for the V0C
     *
     * @ingroup alice_rivet
     */
    const Cut V0Cacceptance = (Cuts::etaIn(-3.7,-1.7)&&(Cuts::abscharge3 > 0));
    /** 
     * The acceptance cut for clusters on layer 0 of the SPD 
     *
     * @ingroup alice_rivet
     */
    const Cut CL0acceptance = (Cuts::etaIn(-2.0,2.0) && (Cuts::abscharge3 > 0));
    /** 
     * The acceptance cut for clusters on layer 1 of the SPD 
     *
     * @ingroup alice_rivet
     */
    const Cut CL1acceptance = (Cuts::etaIn(-1.4,1.4) && (Cuts::abscharge3 > 0));
    /** 
     * The acceptance cut for mid-rapidity 
     *
     * @ingroup alice_rivet
     */
    const Cut Eta1acceptance = (Cuts::etaIn(-1,1) && (Cuts::abscharge3 > 0));
    /** 
     * The acceptance cut for SPD FASTOR 
     *
     * @ingroup alice_rivet
     */
    const Cut FASTORacceptance = CL0acceptance;
#if 0
    /** 
     * Pb identification 
     *
     * @ingroup alice_rivet
     */
    const int PbId = (1000000000 + // ION identifier
		      0*10000000 + // # strange quarks
		        82*10000 + // atomic number
		          208*10 + // atomic weight
		             0*1); // Isomer number 
#endif
  }
}

#endif
//
// EOF
//

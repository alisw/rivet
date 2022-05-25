//FJSTARTHEADER
// $Id: PxConePlugin.cc 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "Rivet/Projections/PxConePlugin.hh"

#include "fastjet/ClusterSequence.hh"
#include <sstream>

// pxcone stuff
// #include "Rivet/Projections/pxcone.h"


// FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace Rivet {

using namespace std;

bool PxConePlugin::_first_time = true;

string PxConePlugin::description () const {
  ostringstream desc;
  
  desc << "PxCone jet algorithm with " 
       << "cone_radius = "        << cone_radius        () << ", "
       << "min_jet_energy = "     << min_jet_energy     () << ", "
       << "overlap_threshold  = " << overlap_threshold  () << ", "
       << "E_scheme_jets  = "     << E_scheme_jets      () 
       << " (NB: non-standard version of PxCone, containing small bug fixes by Gavin Salam)";

  return desc.str();
}


void PxConePlugin::run_clustering(fastjet::ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  //_print_banner(clust_seq.fastjet_banner_stream());
 
  // only have hh mode
  int mode = 2;

  int    ntrak = clust_seq.jets().size(), itkdm = 4;
  double *ptrak = new double[ntrak*4+1];
  for (int i = 0; i < ntrak; i++) {
    ptrak[4*i+0] = clust_seq.jets()[i].px();
    ptrak[4*i+1] = clust_seq.jets()[i].py();
    ptrak[4*i+2] = clust_seq.jets()[i].pz();
    ptrak[4*i+3] = clust_seq.jets()[i].E();
  }  

  // max number of allowed jets
  int mxjet = ntrak;
  int njet;
  double *pjet  = new double[mxjet*5+1];
  int    *ipass = new int[ntrak+1];
  int    *ijmul = new int[mxjet+1];
  int ierr;

  // run pxcone
  pxcone_(
    mode   ,    // 1=>e+e-, 2=>hadron-hadron
    ntrak  ,    // Number of particles
    itkdm  ,    // First dimension of PTRAK array: 
    ptrak  ,    // Array of particle 4-momenta (Px,Py,Pz,E)
    cone_radius()  ,    // Cone size (half angle) in radians
    min_jet_energy() ,    // Minimum Jet energy (GeV)
    overlap_threshold()  ,    // Maximum fraction of overlap energy in a jet
    mxjet  ,    // Maximum possible number of jets
    njet   ,    // Number of jets found
    pjet ,       // 5-vectors of jets
    ipass,      // Particle k belongs to jet number IPASS(k)-1
                // IPASS = -1 if not assosciated to a jet
    ijmul,      // Jet i contains IJMUL[i] particles
    &ierr        // = 0 if all is OK ;   = -1 otherwise
    );

  if (ierr != 0) throw fastjet::Error("An error occurred while running PXCONE");

  // now transfer information back 
  valarray<int> last_index_created(njet);

  vector<vector<int> > jet_particle_content(njet);

  // get a list of particles in each jet
  for (int itrak = 0; itrak < ntrak; itrak++) {
    int jet_i = ipass[itrak] - 1;
    if (jet_i >= 0) jet_particle_content[jet_i].push_back(itrak);
  }

  // now transfer the jets back into our own structure -- we will
  // mimic the cone code with a sequential recombination sequence in
  // which the jets are built up by adding one particle at a time
  for(int ipxjet = njet-1; ipxjet >= 0; ipxjet--) {
    const vector<int> & jet_trak_list = jet_particle_content[ipxjet];
    int jet_k = jet_trak_list[0];
  
    for (unsigned ilist = 1; ilist < jet_trak_list.size(); ilist++) {
      int jet_i = jet_k;
      // retrieve our misappropriated index for the jet
      int jet_j = jet_trak_list[ilist];
      // do a fake recombination step with dij=0
      double dij = 0.0;
      //clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
      if (ilist != jet_trak_list.size()-1 || E_scheme_jets()) {
        // our E-scheme recombination in cases where it doesn't matter
        clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
      } else {
        // put in pxcone's momentum for the last recombination so that the
        // final inclusive jet corresponds exactly to PXCONE's
        clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, 
                      fastjet::PseudoJet(pjet[5*ipxjet+0],pjet[5*ipxjet+1],
                                          pjet[5*ipxjet+2],pjet[5*ipxjet+3]),
                                                 jet_k);
      }
    }
  
    // NB: put a sensible looking d_iB just to be nice...
    double d_iB = clust_seq.jets()[jet_k].perp2();
    clust_seq.plugin_record_iB_recombination(jet_k, d_iB);
  }


  //// following code is for testing only
  //cout << endl;
  //for (int ijet = 0; ijet < njet; ijet++) {
  //  PseudoJet jet(pjet[ijet][0],pjet[ijet][1],pjet[ijet][2],pjet[ijet][3]);
  //  cout << jet.perp() << " " << jet.rap() << endl;
  //}
  //cout << "-----------------------------------------------------\n";
  //vector<PseudoJet> ourjets(clust_seq.inclusive_jets());
  //for (vector<PseudoJet>::const_iterator ourjet = ourjets.begin();
  //     ourjet != ourjets.end(); ourjet++) {
  //  cout << ourjet->perp() << " " << ourjet->rap() << endl;
  //}
  ////cout << endl;

  delete[] ptrak;
  delete[] ipass;
  delete[] ijmul;
  delete[] pjet;
}

// print a banner for reference to the 3rd-party code
void PxConePlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the PxCone plugin for FastJet                           " << endl;
  (*ostr) << "# Original code by the Luis Del Pozo, David Ward and Michael H. Seymour   " << endl;
  (*ostr) << "# If you use this plugin, please cite                                     " << endl;
  (*ostr) << "#   M. H. Seymour and C. Tevlin, JHEP 0611 (2006) 052 [hep-ph/0609100].   " << endl;
  (*ostr) << "# in addition to the usual FastJet reference.                             " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

// FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

/* pxcone.f -- translated by f2c and hacked by Leif Lönnblad to avoid
   linking with libf2c.
*/

//#include "Rivet/Projections/pxcone.h"
using namespace std;

/* Table of constant values, which are actually non const to be able
   to be used as fortran arguments. */
static int MAXV = 20000;
static int VDIM = 3;

void pxtry_(int, double *, int,  double *, double *, double *, double *, 
	    double *, int *, int *);

void pxsorv_(int, double *, int *, char);

void pxsear_(int, double *, int, double *, double *, double *, int &, int *, 
             double *, int *, int *);

void pxolap_(int, int, int, int *, double *, double *, double);

void pxnorv_(int *, double *, double *, int *);

// The standard fortran SIGN function for doubles.
inline double d_sign(double a, double b) {
  return b < 0.0? -fabs(a): fabs(a);
}

// The standard fortran MOD function for doubles.
inline double d_mod(double a, double p) {
  return a - int(a/p)*p;
}

/* ---RETURNS PHI, MOVED ONTO THE RANGE [-PI,PI) */
inline double pxmdpi(double phi) {
  while ( phi <= -M_PI ) phi += 2*M_PI;
  while ( phi > M_PI ) phi -= 2*M_PI;
  return abs(phi) < 1e-15? 0.0: phi;
}
//   if (phi <= M_PI) {
//     if (phi > -M_PI)
//       return abs(phi) < 1e-15? 0.0: phi;
//     else if (phi > -3*M_PI)
//       phi += 2*M_PI;
//     else
//       phi = -d_mod(M_PI - phi, 2*M_PI) + M_PI;
//   } else if (phi <= 3.0*M_PI) {
//     phi -= 2*M_PI;
//   } else {
//     phi = d_mod(phi + M_PI, 2*M_PI) - M_PI;
//   }

//   return abs(phi) < 1e-15? 0.0: phi;

// }

/* Set integer vector a to zero */
inline void pxzeri(int n, int *a){
  for (int i = 0; i < n; ++i) a[i] = 0;
}

/* Set vector a to zero */
inline void pxzerv(int n, double *a) {
    for (int i = 0; i < n; ++i)	a[i] = 0.;
}

/* add vectors c = a + b */
inline void pxaddv(int n, double *a, double *b, double *c) {
  for (int i = 0; i < n; ++i) c[i] = a[i] + b[i];
}

bool pxuvec(int ntrak, double *pp, double *pu) {

  /* Parameter adjustments */
  pu -= 4;
  pp -= 5;

  for (int n = 1; n <= ntrak; ++n) {
    double mag = 0.0;
    for ( int mu = 1; mu <= 3; ++mu)
      mag += pp[mu + (n << 2)]*pp[mu + (n << 2)];
    mag = sqrt(mag);
    if (mag == 0.0 ) {
      printf(" PXCONE: An input particle has zero mod(p)\n");
      return false;
    }
    for (int mu = 1; mu <= 3; ++mu)
      pu[mu + n * 3] = pp[mu + (n << 2)] / mag;
  }
  return true;
}

/* calculate angle between two vectors */
void pxang3(double *a, double *b, double &cost, double &thet) {

  cost = 1.0;
  thet = 0.0;
  double c = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])*
             (b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  if (c <= 0.) return;
  
  c = 1/sqrt(c);
  cost = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*c;
  thet = acos(cost);
  
}

/* ** Note that although JETLIS is assumed to be a 2d array, it */
/* ** it is used as 1d in this routine for efficiency */
/* ** Checks to see if TSTLIS entries correspond to a jet already found */
/* ** and entered in JETLIS */
int pxnew(int *tstlis, int *jetlis, int ntrak, int njet) {

  int match;
  for (int i = 0; i < njet; ++i) {
    match = true;
    int in = i - 5000;
    for (int n = 0; n < ntrak; ++n) {
      in += 5000;
      if (tstlis[n] != jetlis[in]) {
        match = false;
        break;
      }
    }
    if (match) return false;
  }
  return true;
}

/* ** Returns T if the first N elements of LIST1 are the same as the */
/* ** first N elements of LIST2. */
bool pxsame(int *list1, int *list2, int n) {
  for (int i = 0; i < n; ++i)
    if (list1[i] != list2[i])
      return false;
  return true;
}

/* ** Routine to put jets into order and eliminate tose less than EPSLON */
/* ** Puts jets in order of energy: 1 = highest energy etc. */
/* ** Then Eliminate jets with energy below EPSLON */
void pxord(double epslon, int & njet, int ntrak,
	 int *jetlis, double *pj)
{
    /* Local variables */
    static int index[5000];
    static double elist[5000], ptemp[20000]	/* was [4][5000] */;
    static int logtmp[25000000]	/* was [5000][5000] */;



/* ** Copy input arrays. */
    /* Parameter adjustments */
    pj -= 5;
    jetlis -= 5001;

    /* Function Body */
    for (int i = 1; i <= njet; ++i) {
      for (int j = 1; j <= 4; ++j) {
        ptemp[j + (i << 2) - 5] = pj[j + (i << 2)];
      }
      for (int j = 1; j <= ntrak; ++j) {
        logtmp[i + j * 5000 - 5001] = jetlis[i + j * 5000];
      }
    }
    for (int i = 1; i <= njet; ++i) {
      elist[i - 1] = pj[(i << 2) + 4];
    }
    
/* ** Sort the energies... */
    pxsorv_(njet, elist, index, 'I');
/* ** Fill PJ and JETLIS according to sort ( sort is in ascending order!!) */
    for (int i = 1; i <= njet; ++i) {
	for (int j = 1; j <= 4; ++j) {
	    pj[j + (i << 2)] = ptemp[j + (index[njet + 1 - i - 1] << 2) - 5];
	}
	for (int j = 1; j <= ntrak; ++j) {
	    jetlis[i + j * 5000] =
              logtmp[index[njet + 1 - i - 1] + j *  5000 - 5001];
	}
    }
/* * Jets are now in order */
/* ** Now eliminate jets with less than Epsilon energy */
    int nold = njet;
    for (int i = 1; i <= nold; ++i) {
	if (pj[(i << 2) + 4] < epslon) {
	    --njet;
	    pj[(i << 2) + 4] = 0.0;
	}
    }
}

// The main PXCONE function.
void pxcone_(int mode, int ntrak, int itkdm, 
	const double *ptrak, double coner, double epslon, double
	ovlim, int mxjet, int & njet, double *pjet, int *
	ipass, int *ijmul, int *ierr)
{
    /* Initialized data */

    static set<double> rold;
    static set<double> epsold;
    static set<double> ovold;

    /* System generated locals */
    int ptrak_dim1, ptrak_offset, i__1, i__2;
    double d__1, d__2, d__3;

    /* Local variables */
    static double cosr, rsep, ppsq, ptsq, cos2r;
    static int i__, j, n;
    static double vseed[3];
    static int iterr;
    static int n1, n2;
    static double pj[20000]	/* was [4][5000] */, pp[20000]	/* was [4][
	    5000] */;
    static int mu;
    static double pu[15000]	/* was [3][5000] */, cosval;
    static int jetlis[25000000]	/* was [5000][5000] */;
    static int unstbl;
    static double vec1[3], vec2[3];

/* .********************************************************* */
/* . ------ */
/* . PXCONE */
/* . ------ */
/* . */
/* . Code downloaded from the following web page */
/* . */
/* .   http://aliceinfo.cern.ch/alicvs/viewvc/JETAN/pxcone.F?view=markup&pathrev=v4-05-04 */
/* . */
/* . on 17/10/2006 by G. Salam. Permission subsequently granted by Michael */
/* . H. Seymour (on behalf of the PxCone authors) for this code to be */
/* . distributed together with FastJet under the terms of the GNU Public */
/* . License v2 (see the file COPYING in the main FastJet directory). */
/* . */
/* .********** Pre Release Version 26.2.93 */
/* . */
/* . Driver for the Cone  Jet finding algorithm of L.A. del Pozo. */
/* . Based on algorithm from D.E. Soper. */
/* . Finds jets inside cone of half angle CONER with energy > EPSLON. */
/* . Jets which receive more than a fraction OVLIM of their energy from */
/* . overlaps with other jets are excluded. */
/* . Output jets are ordered in energy. */
/* . If MODE.EQ.2 momenta are stored as (eta,phi,<empty>,pt) */
/* . Usage     : */
/* . */
/* .      INTEGER  ITKDM,MXTRK */
/* .      PARAMETER  (ITKDM=4.or.more,MXTRK=1.or.more) */
/* .      INTEGER  MXJET, MXTRAK, MXPROT */
/* .      PARAMETER  (MXJET=10,MXTRAK=500,MXPROT=500) */
/* .      INTEGER  IPASS (MXTRAK),IJMUL (MXJET) */
/* .      INTEGER  NTRAK,NJET,IERR,MODE */
/* .      DOUBLE PRECISION  PTRAK (ITKDM,MXTRK),PJET (5,MXJET) */
/* .      DOUBLE PRECISION  CONER, EPSLON, OVLIM */
/* .      NTRAK = 1.to.MXTRAK */
/* .      CONER   = ... */
/* .      EPSLON  = ... */
/* .      OVLIM   = ... */
/* .      CALL PXCONE (MODE,NTRAK,ITKDM,PTRAK,CONER,EPSLON,OVLIM,MXJET, */
/* .     +             NJET,PJET,IPASS,IJMUL,IERR) */
/* . */
/* . INPUT     :  MODE      1=>e+e-, 2=>hadron-hadron */
/* . INPUT     :  NTRAK     Number of particles */
/* . INPUT     :  ITKDM     First dimension of PTRAK array */
/* . INPUT     :  PTRAK     Array of particle 4-momenta (Px,Py,Pz,E) */
/* . INPUT     :  CONER     Cone size (half angle) in radians */
/* . INPUT     :  EPSLON    Minimum Jet energy (GeV) */
/* . INPUT     :  OVLIM     Maximum fraction of overlap energy in a jet */
/* . INPUT     :  MXJET     Maximum possible number of jets */
/* . OUTPUT    :  NJET      Number of jets found */
/* . OUTPUT    :  PJET      5-vectors of jets */
/* . OUTPUT    :  IPASS(k)  Particle k belongs to jet number IPASS(k) */
/* .                        IPASS = -1 if not assosciated to a jet */
/* . OUTPUT    :  IJMUL(i)  Jet i contains IJMUL(i) particles */
/* . OUTPUT    :  IERR      = 0 if all is OK ;   = -1 otherwise */
/* . */
/* . CALLS     : PXSEAR, PXSAME, PXNEW, PXTRY, PXORD, PXUVEC, PXOLAP */
/* . CALLED    : User */
/* . */
/* . AUTHOR    :  L.A. del Pozo */
/* . CREATED   :  26-Feb-93 */
/* . LAST MOD  :   2-Mar-93 */
/* . */
/* . Modification Log. */
/* . 25-Feb-07: G P Salam   - fix bugs concerning 2pi periodicity in eta phi mode */
/* .                        - added commented code to get consistent behaviour */
/* .                          regardless of particle order (replaces n-way */
/* .                          midpoints with 2-way midpoints however...) */
/* . 2-Jan-97: M Wobisch    - fix bug concerning COS2R in eta phi mode */
/* . 4-Apr-93: M H Seymour  - Change 2d arrays to 1d in PXTRY & PXNEW */
/* . 2-Apr-93: M H Seymour  - Major changes to add boost-invariant mode */
/* . 1-Apr-93: M H Seymour  - Increase all array sizes */
/* . 30-Mar-93: M H Seymour - Change all REAL variables to DOUBLE PRECISION */
/* . 30-Mar-93: M H Seymour - Change OVLIM into an input parameter */
/* . 2-Mar-93: L A del Pozo - Fix Bugs in PXOLAP */
/* . 1-Mar-93: L A del Pozo - Remove Cern library routine calls */
/* . 1-Mar-93: L A del Pozo - Add Print out of welcome and R and Epsilon */
/* . */
/* .********************************************************* */
/* +SEQ,DECLARE. */
/* ** External Arrays */
/* ** Internal Arrays */
/* ** Used in the routine. */
/* MWobisch */
/* MWobisch */
    /* Parameter adjustments */
    --ipass;
    ptrak_dim1 = itkdm;
    ptrak_offset = 1 + ptrak_dim1 * 1;
    ptrak -= ptrak_offset;
    --ijmul;
    pjet -= 6;

    /* Function Body */
/* MWobisch */
/* *************************************** */
    rsep = 2.;
/* *************************************** */
/* MWobisch */
    *ierr = 0;

/* ** INITIALIZE */

/* ** Print welcome and Jetfinder parameters */
    if ((rold.find(coner) == rold.end() ||
         epsold.find(epslon) == epsold.end() ||
         ovold.find(ovlim) == ovold.end()) ) {
      printf("%s\n", " *********** PXCONE: Cone Jet-finder ***********");
      printf("%s\n", "    Written by Luis Del Pozo of OPAL");
      printf("%s\n", "    Modified for eta-phi by Mike Seymour");
      printf("%s\n", "    Includes bug fixes by Wobisch, Salam");
      printf("%s\n", "    Translated to c(++) by Leif Lonnblad");
      printf("%s%5.2f%s\n", "    Cone Size R = ",coner," Radians");
      printf("%s%5.2f%s\n", "    Min Jet energy Epsilon = ",epslon," GeV");
      printf("%s%5.2f\n", "    Overlap fraction parameter = ",ovlim);
      printf("%s\n", "    PXCONE is not a supported product and is");
      printf("%s\n", "    is provided for comparative purposes only");
      printf("%s\n", " ***********************************************");

      rold.insert(coner);
      epsold.insert(epslon);
      ovold.insert(ovlim);
    }

/* ** Copy calling array PTRAK  to internal array PP(4,NTRAK) */

    if (ntrak > 5000) {
/*         WRITE (6,*) ' PXCONE: Ntrak too large: ',NTRAK */
      printf("%s%d\n", " PXCONE: Ntrak too large: ", ntrak);
	*ierr = -1;
	return;
    }
    if (mode != 2) {
	i__1 = ntrak;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		pp[j + (i__ << 2) - 5] = ptrak[j + i__ * ptrak_dim1];
	    }
	}
    } else {
/* ** Converting to eta,phi,pt if necessary */
	i__1 = ntrak;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = ptrak[i__ * ptrak_dim1 + 1];
/* Computing 2nd power */
	    d__2 = ptrak[i__ * ptrak_dim1 + 2];
	    ptsq = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
	    d__3 = ptrak[i__ * ptrak_dim1 + 3];
/* Computing 2nd power */
	    d__2 = sqrt(ptsq + d__3 * d__3) + (d__1 = ptrak[i__ * ptrak_dim1 
		    + 3], abs(d__1));
	    ppsq = d__2 * d__2;
	    if (ptsq <= ppsq * (float)4.25e-18) {
		pp[(i__ << 2) - 4] = 20.;
	    } else {
		pp[(i__ << 2) - 4] = log(ppsq / ptsq) * (float).5;
	    }
	    pp[(i__ << 2) - 4] = d_sign(pp[(i__ << 2) - 4], ptrak[i__ * 
		    ptrak_dim1 + 3]);
	    if (ptsq == 0.) {
		pp[(i__ << 2) - 3] = 0.;
	    } else {
		pp[(i__ << 2) - 3] = atan2(ptrak[i__ * ptrak_dim1 + 2], ptrak[
			i__ * ptrak_dim1 + 1]);
	    }
	    pp[(i__ << 2) - 2] = 0.;
	    pp[(i__ << 2) - 1] = sqrt(ptsq);
	    pu[i__ * 3 - 3] = pp[(i__ << 2) - 4];
	    pu[i__ * 3 - 2] = pp[(i__ << 2) - 3];
	    pu[i__ * 3 - 1] = pp[(i__ << 2) - 2];
	}
    }

/* ** Zero output variables */

    njet = 0;
    i__1 = ntrak;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 5000; ++j) {
	    jetlis[j + i__ * 5000 - 5001] = false;
	}
    }
    pxzerv(MAXV, pj);
    pxzeri(mxjet, &ijmul[1]);

    if (mode != 2) {
	cosr = cos(coner);
	cos2r = cos(coner);
    } else {
/* ** Purely for convenience, work in terms of 1-R**2 */
/* Computing 2nd power */
	d__1 = coner;
	cosr = 1 - d__1 * d__1;
/* MW -- select Rsep: 1-(Rsep*CONER)**2 */
/* Computing 2nd power */
	d__1 = rsep * coner;
	cos2r = 1 - d__1 * d__1;
/* ORIGINAL         COS2R =  1-(2*CONER)**2 */
    }
    unstbl = false;
    if (mode != 2) {
      if ( !pxuvec(ntrak, pp, pu) ) {
        *ierr = 1;
        return;
      }
    }
/* ** Look for jets using particle diretions as seed axes */

    i__1 = ntrak;
    for (n = 1; n <= i__1; ++n) {
	for (mu = 1; mu <= 3; ++mu) {
	    vseed[mu - 1] = pu[mu + n * 3 - 4];
	}
	pxsear_(mode, &cosr, ntrak, pu, pp, vseed, njet, jetlis, pj, &unstbl, 
		ierr);
	if (*ierr != 0) {
	    return;
	}
    }

    i__1 = njet - 1;
    for (n1 = 1; n1 <= i__1; ++n1) {
	vec1[0] = pj[(n1 << 2) - 4];
	vec1[1] = pj[(n1 << 2) - 3];
	vec1[2] = pj[(n1 << 2) - 2];
	if (mode != 2) {
	    pxnorv_(&VDIM, vec1, vec1, &iterr);
	}
/*         DO 150 N2 = N1+1,NJTORG ! GPS -- to get consistent behaviour */
	i__2 = njet;
	for (n2 = n1 + 1; n2 <= i__2; ++n2) {
	    vec2[0] = pj[(n2 << 2) - 4];
	    vec2[1] = pj[(n2 << 2) - 3];
	    vec2[2] = pj[(n2 << 2) - 2];
	    if (mode != 2) {
		pxnorv_(&VDIM, vec2, vec2, &iterr);
	    }
	    pxaddv(VDIM, vec1, vec2, vseed);
	    if (mode != 2) {
		pxnorv_(&VDIM, vseed, vseed, &iterr);
	    } else {
		vseed[0] /= 2;
/* VSEED(2)=VSEED(2)/2 */
/* GPS 25/02/07 */
		d__2 = vec2[1] - vec1[1];
		d__1 = vec1[1] + pxmdpi(d__2) * .5;
		vseed[1] = pxmdpi(d__1);
	    }
/* ---ONLY BOTHER IF THEY ARE BETWEEN 1 AND 2 CONE RADII APART */
	    if (mode != 2) {
		cosval = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * 
			vec2[2];
	    } else {
		if (abs(vec1[0]) >= 20. || abs(vec2[0]) >= 20.) {
		    cosval = -1e3;
		} else {
/* Computing 2nd power */
		    d__1 = vec1[0] - vec2[0];
		    d__3 = vec1[1] - vec2[1];
/* Computing 2nd power */
		    d__2 = pxmdpi(d__3);
		    cosval = 1 - (d__1 * d__1 + d__2 * d__2);
		}
	    }
	    if (cosval <= cosr && cosval >= cos2r) {
		pxsear_(mode, &cosr, ntrak, pu, pp, vseed, njet, jetlis, pj, &
			unstbl, ierr);
	    }
/*            CALL PXSEAR(MODE,COSR,NTRAK,PU,PP,VSEED,NJET, */
/*     +           JETLIS,PJ,UNSTBL,IERR) */
	    if (*ierr != 0) {
		return;
	    }
	}
    }
    if (unstbl) {
	*ierr = -1;
/*        WRITE (6,*) ' PXCONE: Too many iterations to find a proto-jet' */
        printf(" PXCONE: Too many iterations to find a proto-jet\n");
	return;
    }
/* ** Now put the jet list into order by jet energy, eliminating jets */
/* ** with energy less than EPSLON. */
    pxord(epslon, njet, ntrak, jetlis, pj);

/* ** Take care of jet overlaps */
    pxolap_(mode, njet, ntrak, jetlis, pj, pp, ovlim);

/* ** Order jets again as some have been eliminated, or lost energy. */
    pxord(epslon, njet, ntrak, jetlis, pj);

/* ** All done!, Copy output into output arrays */
    if (njet > mxjet) {
/*         WRITE (6,*) ' PXCONE:  Found more than MXJET jets' */
      printf(" PXCONE:  Found more than MXJET jets\n");
	*ierr = -1;
	goto L99;
    }
    if (mode != 2) {
	i__1 = njet;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		pjet[j + i__ * 5] = pj[j + (i__ << 2) - 5];
	    }
	}
    } else {
	i__1 = njet;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pjet[i__ * 5 + 1] = pj[(i__ << 2) - 1] * cos(pj[(i__ << 2) - 3]);
	    pjet[i__ * 5 + 2] = pj[(i__ << 2) - 1] * sin(pj[(i__ << 2) - 3]);
	    pjet[i__ * 5 + 3] = pj[(i__ << 2) - 1] * sinh(pj[(i__ << 2) - 4]);
	    pjet[i__ * 5 + 4] = pj[(i__ << 2) - 1] * cosh(pj[(i__ << 2) - 4]);
	}
    }
    i__1 = ntrak;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipass[i__] = -1;
	i__2 = njet;
	for (j = 1; j <= i__2; ++j) {
	    if (jetlis[j + i__ * 5000 - 5001]) {
		++ijmul[j];
		ipass[i__] = j;
	    }
	}
    }
L99:
    return;
} /* pxcone_ */


void pxnorv_(int *n, double *a, double *b, int *iterr)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */

    /* Local variables */
    static double c__;
    static int i__;

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    c__ = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = a[i__];
	c__ += d__1 * d__1;
    }
    if (c__ <= 0.) {
	return;
    }
    c__ = 1 / sqrt(c__);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = a[i__] * c__;
    }
    return;
} /* pxnorv_ */


void pxolap_(int mode, int njet, int ntrak, 
	int *jetlis, double *pj, double *pp, double ovlim)
{
    /* Initialized data */

    static int ijmin = 0;

    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1, d__2, d__3;

    /* Local variables */
    static int ijet[5000];
    static double thet, cost;
    static int i__, j, n;
    static double eover, thmin;
    static int nj, mu;
    static int ovelap;
    static double vec1[3], vec2[3];


/* ** Looks for particles assigned to more than 1 jet, and reassigns them */
/* ** If more than a fraction OVLIM of a jet's energy is contained in */
/* ** higher energy jets, that jet is neglected. */
/* ** Particles assigned to the jet closest in angle (a la CDF, Snowmass). */
/* +SEQ,DECLARE. */
    /* Parameter adjustments */
    pp -= 5;
    pj -= 5;
    jetlis -= 5001;

    /* Function Body */

    if (njet <= 1) {
	return;
    }
/* ** Look for jets with large overlaps with higher energy jets. */
    i__1 = njet;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* ** Find overlap energy between jets I and all higher energy jets. */
	eover = (float)0.;
	i__2 = ntrak;
	for (n = 1; n <= i__2; ++n) {
	    ovelap = false;
	    i__3 = i__ - 1;
	    for (j = 1; j <= i__3; ++j) {
		if (jetlis[i__ + n * 5000] && jetlis[j + n * 5000]) {
		    ovelap = true;
		}
	    }
	    if (ovelap) {
		eover += pp[(n << 2) + 4];
	    }
	}
/* ** Is the fraction of energy shared larger than OVLIM? */
	if (eover > ovlim * pj[(i__ << 2) + 4]) {
/* ** De-assign all particles from Jet I */
	    i__2 = ntrak;
	    for (n = 1; n <= i__2; ++n) {
		jetlis[i__ + n * 5000] = false;
	    }
	}
    }
/* ** Now there are no big overlaps, assign every particle in */
/* ** more than 1 jet to the closet jet. */
/* ** Any particles now in more than 1 jet are assigned to the CLOSET */
/* ** jet (in angle). */
    i__1 = ntrak;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nj = 0;
	i__2 = njet;
	for (j = 1; j <= i__2; ++j) {
	    if (jetlis[j + i__ * 5000]) {
		++nj;
		ijet[nj - 1] = j;
	    }
	}
	if (nj > 1) {
/* ** Particle in > 1 jet - calc angles... */
	    vec1[0] = pp[(i__ << 2) + 1];
	    vec1[1] = pp[(i__ << 2) + 2];
	    vec1[2] = pp[(i__ << 2) + 3];
	    thmin = (float)0.;
	    i__2 = nj;
	    for (j = 1; j <= i__2; ++j) {
		vec2[0] = pj[(ijet[j - 1] << 2) + 1];
		vec2[1] = pj[(ijet[j - 1] << 2) + 2];
		vec2[2] = pj[(ijet[j - 1] << 2) + 3];
		if (mode != 2) {
		    pxang3(vec1, vec2, cost, thet);
		} else {
/* Computing 2nd power */
		    d__1 = vec1[0] - vec2[0];
		    d__3 = vec1[1] - vec2[1];
/* Computing 2nd power */
		    d__2 = pxmdpi(d__3);
		    thet = d__1 * d__1 + d__2 * d__2;
		}
		if (j == 1) {
		    thmin = thet;
		    ijmin = ijet[j - 1];
		} else if (thet < thmin) {
		    thmin = thet;
		    ijmin = ijet[j - 1];
		}
	    }
/* ** Assign track to IJMIN */
	    i__2 = njet;
	    for (j = 1; j <= i__2; ++j) {
		jetlis[j + i__ * 5000] = false;
	    }
	    jetlis[ijmin + i__ * 5000] = true;
	}
    }
/* ** Recompute PJ */
    i__1 = njet;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (mu = 1; mu <= 4; ++mu) {
	    pj[mu + (i__ << 2)] = (float)0.;
	}
	i__2 = ntrak;
	for (n = 1; n <= i__2; ++n) {
	    if (jetlis[i__ + n * 5000]) {
		if (mode != 2) {
		    for (mu = 1; mu <= 4; ++mu) {
			pj[mu + (i__ << 2)] += pp[mu + (n << 2)];
		    }
		} else {
		    pj[(i__ << 2) + 1] += pp[(n << 2) + 4] / (pp[(n << 2) + 4]
			     + pj[(i__ << 2) + 4]) * (pp[(n << 2) + 1] - pj[(
			    i__ << 2) + 1]);
/* GPS 25/02/07 */
		    d__2 = pp[(n << 2) + 2] - pj[(i__ << 2) + 2];
		    d__1 = pj[(i__ << 2) + 2] + pp[(n << 2) + 4] / (pp[(n << 
			    2) + 4] + pj[(i__ << 2) + 4]) * pxmdpi(d__2);
		    pj[(i__ << 2) + 2] = pxmdpi(d__1);
/*                PJ(2,I)=PJ(2,I) */
/*     +               + PP(4,N)/(PP(4,N)+PJ(4,I))*PXMDPI(PP(2,N)-PJ(2,I)) */
		    pj[(i__ << 2) + 4] += pp[(n << 2) + 4];
		}
	    }
	}
    }
    return;
} /* pxolap_ */


/* ******************************************************************* */
void pxsear_(int mode, double *cosr, int ntrak, 
             double *pu, double *pp, double *vseed, int & njet, 
             int *jetlis, double *pj, int *unstbl, int *ierr)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int iter;
    static double pnew[4];
    static int n;
    static double naxis[3], oaxis[3];
    static int ok;
    static int mu;
    static int oldlis[5000];
    static int newlis[5000];


/* +SEQ,DECLARE. */
/* ** Using VSEED as a trial axis , look for a stable jet. */
/* ** Check stable jets against those already found and add to PJ. */
/* ** Will try up to MXITER iterations to get a stable set of particles */
/* ** in the cone. */

    /* Parameter adjustments */
    pj -= 5;
    jetlis -= 5001;
    --vseed;
    pp -= 5;
    pu -= 4;

    /* Function Body */
    for (mu = 1; mu <= 3; ++mu) {
	oaxis[mu - 1] = vseed[mu];
    }
    i__1 = ntrak;
    for (n = 1; n <= i__1; ++n) {
	oldlis[n - 1] = false;
    }
    for (iter = 1; iter <= 30; ++iter) {
	pxtry_(mode, cosr, ntrak, &pu[4], &pp[5], oaxis, naxis, pnew, newlis, 
		&ok);
/* ** Return immediately if there were no particles in the cone. */
	if (! ok) {
	    return;
	}
	if (pxsame(newlis, oldlis, ntrak)) {
/* ** We have a stable jet. */
	    if (pxnew(newlis, &jetlis[5001], ntrak, njet)) {
/* ** And the jet is a new one. So add it to our arrays. */
/* ** Check arrays are big anough... */
		if (njet == 5000) {
/*             WRITE (6,*) ' PXCONE:  Found more than MXPROT proto-jets' */
                  printf(" PXCONE:  Found more than MXPROT proto-jets\n");
		    *ierr = -1;
		    return;
		}
		++njet;
		i__1 = ntrak;
		for (n = 1; n <= i__1; ++n) {
		    jetlis[njet + n * 5000] = newlis[n - 1];
		}
		for (mu = 1; mu <= 4; ++mu) {
		    pj[mu + (njet << 2)] = pnew[mu - 1];
		}
	    }
	    return;
	}
/* ** The jet was not stable, so we iterate again */
	i__1 = ntrak;
	for (n = 1; n <= i__1; ++n) {
	    oldlis[n - 1] = newlis[n - 1];
	}
	for (mu = 1; mu <= 3; ++mu) {
	    oaxis[mu - 1] = naxis[mu - 1];
	}
    }
    *unstbl = true;
    return;
} /* pxsear_ */


void pxsorv_(int n, double *a, int *k, char opt)
{
    /* System generated locals */
    int i__1;


    /* Local variables */
    static double b[5000];
    static int i__, j, il[5000], ir[5000];

/*     Sort A(N) into ascending order */
/*     OPT = 'I' : return index array K only */
/*     OTHERWISE : return sorted A and index array K */
/* ----------------------------------------------------------------------- */

/*      INT N,I,J,K(N),IL(NMAX),IR(NMAX) */
/* LUND */

/*      DOUBLE PRECISION A(N),B(NMAX) */
/* LUND */
    /* Parameter adjustments */
    --k;
    --a;

    /* Function Body */
    if (n > 5000) {
      // WRITE	s_stop("Sorry, not enough room in Mike's PXSORV", (ftnlen)39);
      printf("Sorry, not enough room in Mike's PXSORV\n");
      abort();
    }
    il[0] = 0;
    ir[0] = 0;
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	il[i__ - 1] = 0;
	ir[i__ - 1] = 0;
	j = 1;
L2:
	if (a[i__] > a[j]) {
	    goto L5;
	}
	if (il[j - 1] == 0) {
	    goto L4;
	}
	j = il[j - 1];
	goto L2;
L4:
	ir[i__ - 1] = -j;
	il[j - 1] = i__;
	goto L10;
L5:
	if (ir[j - 1] <= 0) {
	    goto L6;
	}
	j = ir[j - 1];
	goto L2;
L6:
	ir[i__ - 1] = ir[j - 1];
	ir[j - 1] = i__;
L10:
	;
    }
    i__ = 1;
    j = 1;
    goto L8;
L20:
    j = il[j - 1];
L8:
    if (il[j - 1] > 0) {
	goto L20;
    }
L9:
    k[i__] = j;
    b[i__ - 1] = a[j];
    ++i__;
    if ((i__1 = ir[j - 1]) < 0) {
	goto L12;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L13;
    }
L13:
    j = ir[j - 1];
    goto L8;
L12:
    j = -ir[j - 1];
    goto L9;
L30:
    if ( opt == 'I') {
	return;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = b[i__ - 1];
    }
    return;
} /* pxsorv_ */

/* ******************************************************************** */

void pxtry_(int mode, double *cosr, int ntrak, 
	double *pu, double *pp, double *oaxis, double *naxis, 
	double *pnew, int *newlis, int *ok)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3;

    /* Builtin functions */

    /* Local variables */
    static double norm;
    static int n, mu;
    static double cosval;
    static double normsq;
    static int npp, npu;


/* +SEQ,DECLARE. */
/* ** Note that although PU and PP are assumed to be 2d arrays, they */
/* ** are used as 1d in this routine for efficiency */
/* ** Finds all particles in cone of size COSR about OAXIS direction. */
/* ** Calculates 4-momentum sum of all particles in cone (PNEW) , and */
/* ** returns this as new jet axis NAXIS (Both unit Vectors) */

    /* Parameter adjustments */
    --newlis;
    --pnew;
    --naxis;
    --oaxis;
    --pp;
    --pu;

    /* Function Body */
    *ok = false;
    for (mu = 1; mu <= 4; ++mu) {
	pnew[mu] = (float)0.;
    }
    npu = -3;
    npp = -4;
    i__1 = ntrak;
    for (n = 1; n <= i__1; ++n) {
	npu += 3;
	npp += 4;
	if (mode != 2) {
	    cosval = (float)0.;
	    for (mu = 1; mu <= 3; ++mu) {
		cosval += oaxis[mu] * pu[mu + npu];
	    }
	} else {
	    if ((d__1 = pu[npu + 1], abs(d__1)) >= 20. || abs(oaxis[1]) >= 
		    20.) {
		cosval = -1e3;
	    } else {
/* Computing 2nd power */
		d__1 = oaxis[1] - pu[npu + 1];
		d__3 = oaxis[2] - pu[npu + 2];
/* Computing 2nd power */
		d__2 = pxmdpi(d__3);
		cosval = 1 - (d__1 * d__1 + d__2 * d__2);
	    }
	}
	if (cosval >= *cosr) {
	    newlis[n] = true;
	    *ok = true;
	    if (mode != 2) {
		for (mu = 1; mu <= 4; ++mu) {
		    pnew[mu] += pp[mu + npp];
		}
	    } else {
		pnew[1] += pp[npp + 4] / (pp[npp + 4] + pnew[4]) * (pp[npp + 
			1] - pnew[1]);
/*                PNEW(2)=PNEW(2) */
/*     +              + PP(4+NPP)/(PP(4+NPP)+PNEW(4)) */
/*     +               *PXMDPI(PP(2+NPP)-PNEW(2)) */
/* GPS 25/02/07 */
		d__2 = pp[npp + 2] - pnew[2];
		d__1 = pnew[2] + pp[npp + 4] / (pp[npp + 4] + pnew[4]) * 
			pxmdpi(d__2);
		pnew[2] = pxmdpi(d__1);
		pnew[4] += pp[npp + 4];
	    }
	} else {
	    newlis[n] = false;
	}
    }
/* ** If there are particles in the cone, calc new jet axis */
    if (*ok) {
	if (mode != 2) {
	    normsq = (float)0.;
	    for (mu = 1; mu <= 3; ++mu) {
/* Computing 2nd power */
		d__1 = pnew[mu];
		normsq += d__1 * d__1;
	    }
	    norm = sqrt(normsq);
	} else {
	    norm = 1.;
	}
	for (mu = 1; mu <= 3; ++mu) {
	    naxis[mu] = pnew[mu] / norm;
	}
    }
} /* pxtry_ */


}

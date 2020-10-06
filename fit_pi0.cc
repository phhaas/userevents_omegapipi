// Code necessary for fitting \pi^0s reconstructed from two clusters
// in ECAL2
// This code is documented in COMPASS note 2007-10
#include "G3part.h"
#include "PaEvent.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#include "kinfit.h"
#include "fit_pi0.h"

/* Coordinates are (x1, y1, E1, x2, y2, E2).  The constraint
   m_\gamma = 0 is enforced by deriving the pz' from these quantities.
   x = X/R, y = Y/R, z = Z/R = sqrt(1-x^2-y^2).  */

/* Squared-mass difference from the \pi0-mass.  */
static void
massdiff (const TVectorD& x, TVectorD& values)
{
  const double& x1 = x[0];
  const double& y1 = x[1];
  const double& E1 = x[2];
  const double& x2 = x[3];
  const double& y2 = x[4];
  const double& E2 = x[5];

  const double z12 = 1 - x1*x1 - y1*y1;
  const double z1 = sqrt (z12);
  const double z22 = 1 - x2*x2 - y2*y2;
  const double z2 = sqrt (z22);

  const double d =
    2*E1*E2*(1 - x1*x2 - y1*y2 - z1*z2) - G3partMass[7]*G3partMass[7];
                                         /* \pi^0 mass-squared.  */

  values.ResizeTo(1);
  values[0] = d;
}


/* Derivative of massdiff.  */

static void
dmassdiff (const TVectorD& x, TMatrixD& values)
{
  const double& x1 = x[0];
  const double& y1 = x[1];
  const double& E1 = x[2];
  const double& x2 = x[3];
  const double& y2 = x[4];
  const double& E2 = x[5];

  const double z12 = 1 - x1*x1 - y1*y1;
  const double z1 = sqrt (z12);
  const double z22 = 1 - x2*x2 - y2*y2;
  const double z2 = sqrt (z22);

  const double d[] = 
    {
      2*E1*E2*(-x2 + z2*x1/z1),
      2*E1*E2*(-y2 + z2*y1/z1),
      2*E2*(1 - x1*x2 - y1*y2 - z1*z2),
      2*E1*E2*(-x1 + z1*x2/z2),
      2*E1*E2*(-y1 + z1*y2/z2),
      2*E1*(1 - x1*x2 - y1*y2 - z1*z2), };

  values.ResizeTo(1, 6);
  values.SetMatrixArray(d);
}


/* Generate covariance matrix for calorimeter cluster, with
   direction and energies of the incoming particle as variables.
   deltaE allows choosing between different error esxtimates for the
   calorimeter clusters.  See inside for specifics.  */
static TMatrixDSym
cov_mat (const PaCaloClus& cl, const PaVertex& v, int deltaE = 0)
{
  /* Linearized change of variables (X, Y, Z, E) -> (x, y, E).
     x,y are the direction cosines, i.e. x = X / R, y = Y / R  */
  double x = cl.X() - v.Pos(0);
  double y = cl.Y() - v.Pos(1);
  double z = cl.Z() - v.Pos(2);
  double x2 = x*x, y2 = y*y, z2 = z*z;
  double r = hypot (hypot (x, y), z);
  double r3 = 1/(r*r*r);

  double a[] =
    { (y2+z2)*r3,    -x*y*r3,    -x*z*r3, 0,
         -x*y*r3, (x2+z2)*r3,    -y*z*r3, 0,
               0,          0,          0, 1, };
  TMatrixD A (3,4,a);

  /* Construct the covariance matrix of (X, Y, Z, E) from the data
     phast gives for the cluster and the vertex.  */

  /* Decide which error estimate to use for the cluster energy.  */
  double sigE2;
  switch(deltaE)
    {
    case 0:
      sigE2 = cl.Eerr() * cl.Eerr();
      break;
    case 1:
      sigE2 = 0.055*0.055*cl.E() + 0.015*0.015*cl.E()*cl.E();
      break;
    case 2:
      sigE2 = 0.07*0.07*cl.E() + 0.015*0.015*cl.E()*cl.E();
      break;
    default:
      abort();
    }

  /* Put together the covariance matrix.  We don't take the
     uncertainty in the original vertex position into account.  This
     affects the directions, but since there are other particles in
     the interaction, a consistent approach would have to refit their
     tracks, which is unwieldy.  Also verified that introducing the
     vertex' uncertainties doesn't significantly alter the
     results.  */
  const float* d = cl.Cov();
  double c[16] =
    { d[0], d[1], d[3], 0,
      d[1], d[2], d[4], 0,
      d[3], d[4], d[5], 0,
      0, 0, 0, sigE2, };
  TMatrixDSym C (4, c);

  // Finally, calculate the covariance matrix for (x, y, E).
  C.Similarity (A);
  return C;
}


// Kinematic fit of the 4-momentum pcorr of a \pi^0
//
// returns true on success, false otherwise
// Input arguments are the two calorimeter clusters, the vertex
// where the \pi^0 is produced, and the optional deltaE which chooses
// between different error estimates for the energies of the
// calorimeter clusters.  See cov_mat() for details.
// Output arguments are the corrected 4-momentum pcorr, the \chi^2 of
// the fit, and the vector of pulls
bool
fit_pi0 (const PaCaloClus& cluster1, const PaCaloClus& cluster2,
	 const PaVertex& vertex, TLorentzVector& pcorr, double& chi2,
	 TVectorD& pulls, int deltaE)
{
  TVector3 v(vertex.Pos(0), vertex.Pos(1), vertex.Pos(2));
  TVector3 v1(cluster1.X() - vertex.Pos(0),
	      cluster1.Y() - vertex.Pos(1),
	      cluster1.Z() - vertex.Pos(2));
  v1 = v1.Unit();
  TVector3 v2(cluster2.X() - vertex.Pos(0),
	      cluster2.Y() - vertex.Pos(1),
	      cluster2.Z() - vertex.Pos(2));
  v2 = v2.Unit();

  TMatrixDSym Cm (6);
  TMatrixDSub (Cm, 0, 2, 0, 2) = cov_mat (cluster1, vertex, deltaE);
  TMatrixDSub (Cm, 3, 5, 3, 5) = cov_mat (cluster2, vertex, deltaE);

  double w[] =
    { v1.X(), v1.Y(), cluster1.E(),
      v2.X(), v2.Y(), cluster2.E(), };
  TVectorD V (6, w);
  bool success = kinfit (V, Cm, massdiff, dmassdiff, pulls);
  TMatrixDSym dcov;
  if (success)
    kinfit_dcov (V, Cm, dmassdiff, dcov);
  else
    return false;
  double newE1 = cluster1.E() + pulls[2];
  double newE2 = cluster2.E() + pulls[5];
  double newpx1 = (v1.X() + pulls[0])*newE1;
  double newpy1 = (v1.Y() + pulls[1])*newE1;
  double newpx2 = (v2.X() + pulls[3])*newE2;
  double newpy2 = (v2.Y() + pulls[4])*newE2; 
  double newpz1 = sqrt(newE1*newE1 - newpx1*newpx1 - newpy1*newpy1);
  double newpz2 = sqrt(newE2*newE2 - newpx2*newpx2 - newpy2*newpy2);
  pcorr = TLorentzVector(newpx1 + newpx2, newpy1 + newpy2,
			 newpz1 + newpz2, newE1 + newE2);

  // Calculate \chi^2
  Cm.Invert();
  // Old ROOT versions don't allow this:
  //  td.chi2 = Cm.Similarity (pulls);
  // hence the following:
  chi2 = 0;
  for (int i=0; i<6; i++)
    for (int j=0; j<6; j++)
      chi2 += pulls[i] * Cm[i][j] * pulls[j];

  // Calculate pull quantities
  for (int i=0; i<6; i++)
    pulls[i] = pulls[i] / sqrt(dcov[i][i]);

  return true;
}

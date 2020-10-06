/* Constrained fit of the measured values X with the covariance matrix
   COV.  The constraints are supposed to be fulfilled where the
   constraint function pF is zero, its derivative is pdF.  The
   iteration is controlled by the parameters precision_goal (stops if
   all components of the step vector are smaller than precision_goal)
   and max_iterations (stops after max_iterations steps).

   Author: Tobias Schlüter  <tobias.schlueter@physik.uni-muenchen.de>

   This code is documented in COMPASS note 2007-10  */

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"

#include "kinfit.h"

static void
kinfit_step (const TVectorD& x, const TMatrixDSym& cov,
	     F_t pF, dF_t pdF, TVectorD& dx)
{
  TVectorD c;
  pF(x, c);
  TMatrixD B;
  pdF(x, B);

  TMatrixDSym BcovBT (cov);
  BcovBT.Similarity (B);

  /* Special casing a one-dimensional constraint (i.e. BcovBT a 1x1
     matrix) brings no significant benefit.  */
  TMatrixDSym BcovBTinv (TDecompChol(BcovBT).Invert());
  TMatrixD A(cov, TMatrixD::kMult,
	     TMatrixD (B, TMatrixD::kTransposeMult,
		       BcovBTinv));
  c *= A;   // c -> A*c
  dx.ResizeTo (c);
  dx = c;
}

static TVectorD x_2nd_to_last;

bool
kinfit (const TVectorD& x0, const TMatrixDSym& cov, F_t F, dF_t dF,
	TVectorD& dx, Double_t precision_goal, Int_t max_iterations)
{
  dx.ResizeTo(x0);
  dx = 0.;

  for (int i=0; i < max_iterations; i++)
    {
      TVectorD b;
      kinfit_step (x0 + dx, cov, F, dF, b);
      if (b.NormInf () < precision_goal)
	{
	  x_2nd_to_last.ResizeTo (x0);
	  x_2nd_to_last = x0 + dx;
	  dx -= b;
	  //cout << "Converged at step " << i+1 << endl;
	  return true;
	}
      dx -= b;
    }

  // didn't converge
  return false;
}


/* New covariance matrix after the constrained fit.  */
void
kinfit_dcov (const TVectorD& x0, const TMatrixDSym& cov, dF_t pdF,
	     TMatrixDSym& dcov)
{
  // Note that the argument x0 is unused.  Instead of using
  // x_2nd_to_last, one could also use x0.  The numerical difference
  // is usually negligible, as the iteration usually converges in very
  // few steps.  In order to not change the calling
  // convention when trying this, the argument was kept.

  // If the user calls this after kinfit returned false, this will
  // return something random, and the user will get what he deserves
  // for not checking return values.

  TMatrixD B;
  pdF(x_2nd_to_last, B);

  TMatrixDSym BcovBT(cov);
  BcovBT.Similarity(B);
  TMatrixDSym dc(TDecompChol(BcovBT).Invert());
  dc.SimilarityT(B);
  dc.Similarity(cov);
  dcov.ResizeTo(dc);
  dcov = dc;
}


#if defined(__GNUC__) && defined (__MAKECINT__) // i.e. if compiling from ROOT
// Due to shortcomings in Cint (as of version 5.16.21), this example only
// works in compiled mode, load with .L kinfit.cc+

#include <math.h>
#include "TH1F.h"
#include "TRandom.h"
#include "TCanvas.h"

///////////////////////////////////////////////////////////////////////////
// Example.

/* Example with two constraints
   Inside angles of a quadrangle

   Let the angles
    ABC, CAB, ABD, BDC, BCD, ACD
   be the measured quantities (in degrees)

   Then the following relations should hold:
    ABC + CAB + ACD - BCD = 180
    ABD - ABC + BDC + BCD = 180  */

static void
H (const TVectorD& x, TVectorD& values)
{
  double d[] =
    { x[0] + x[1] + x[5] - x[4] - 180,
      x[2] - x[0] + x[3] + x[4] - 180 };

  values.ResizeTo(2);
  values[0] = d[0];
  values[1] = d[1];
}


static void
dH (const TVectorD& x, TMatrixD& values)
{
  double m[] =
    {  1, 1, 0, 0, -1, 1,
       -1, 0, 1, 1, 1, 0, };

  values.ResizeTo(2, 6);
  values.SetMatrixArray(m);
}

void
example (int n = 10000, double error_wrong_by = 0)
{
  TH1F* hconf = new TH1F("hconf", "confidence level distribution",
			 50, 0, 1);
  TH1F* hpull = new TH1F("hpull", "pull example", 100, -4, 4);

  // Setup the covariance matrix and its inverse, the weight matrix
  double m[] = { 9, 0, 0, 0, 0, 0,
		 0, 4, 0, 0, 0, 0,
		 0, 0, 5, 0, 0, 0,
		 0, 0, 0, 4, 0, 0,
		 0, 0, 0, 0, 4, 0,
                 0, 0, 0, 0, 0, 4 };
  TMatrixDSym Cm(6, m);
  TMatrixDSym Gm(Cm);
  Gm.Invert ();

  // Generate a series of fits
  for (int i=0; i<n; i++)
    {
      TRandom* r = gRandom;
      double v[] = { r->Gaus(70,3), r->Gaus(100,2),
		     r->Gaus(110,sqrt(5)), r->Gaus(130,2),
		     r->Gaus(10,2), r->Gaus(20,2+error_wrong_by) };

      TVectorD x (6, v);

      TVectorD dx;
      kinfit (x, Cm, H, dH, dx);
      TMatrixDSym dCm;
      kinfit_dcov(x, Cm, dH, dCm);
      hpull->Fill(dx[5] / sqrt(dCm[5][5]));
      // Replace the following with \sum_i\sum_j dx_i Gm_{ij} dx_j in
      // ROOT versions < 5
      double chi2 = Gm.Similarity(dx);
      hconf->Fill(exp(-chi2/2));
    }

  TCanvas *c = new TCanvas ("c", "c");
  c->Divide(2);
  c->cd(1);
  hconf->Draw();
  hconf->Fit("pol1");
  c->cd(2);
  hpull->Draw();
  hpull->Fit("gaus");
}


#endif // __CINT__

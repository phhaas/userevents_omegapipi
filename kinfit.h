/* This code is documented in COMPASS note 2007-10  */
#ifndef KINFIT_H
#define KINFIT_H

//#include "TVectorD.h"
//#include "TMatrixD.h"

typedef void (* F_t) (const TVectorD& x, TVectorD& values);
typedef void (*dF_t) (const TVectorD& x, TMatrixD& values);

bool kinfit (const TVectorD& x0, const TMatrixDSym& cov,
	     F_t F, dF_t dF, TVectorD& dx,
	     Double_t precision_goal = 1e-12, Int_t max_iterations = 10);
void kinfit_dcov (const TVectorD& x0, const TMatrixDSym& cov,
		  dF_t pdF, TMatrixDSym& dcov);

#endif

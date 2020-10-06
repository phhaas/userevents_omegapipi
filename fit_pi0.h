// This code is documented in COMPASS note 2007-10
#ifndef FIT_MASS
#define FIT_MASS

#include "PaEvent.h"
#include "TLorentzVector.h"
#include "TVector.h"

bool
fit_pi0 (const PaCaloClus& cluster1, const PaCaloClus& cluster2,
	 const PaVertex& vertex, TLorentzVector& pcorr, double& chi2,
	 TVectorD& pulls, int deltaE = 0);

#endif

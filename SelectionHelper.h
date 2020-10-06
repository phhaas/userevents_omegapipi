#ifndef SELECTIONHELPER_H
#define SELECTIONHELPER_H

#include "TH1D.h"

#include "Phast.h"

using namespace std;

class SelectionHelper
{

public:

	SelectionHelper(TTree& tree);

	virtual ~SelectionHelper() { }

	void reset() { setIsAcceptedMC(false); }
	void fill           (const PaEvent& event);
	void fillStatistics (const TString& name) { _hStatistics->Fill(name, 1); }
	void handleMCEvent  (TTree& tree, const bool accepted, const PaEvent& event);
	bool isCorrectCharge(const PaParticle& beamParticle, const vector<PaParticle>& particles) const;
	bool isDT0() const;

	static bool isInRange(const double value,
	                      const double lowerBound,
	                      const double upperBound)
	{ return lowerBound <= value and value <= upperBound; }

private:

	void setIsAcceptedMC(const int  accepted) { _isAcceptedMC = accepted;       }
	void setIsAcceptedMC(const bool accepted) { setIsAcceptedMC((int)accepted); }

	bool hasTrackParameters(const PaParticle& particle) const { return particle.NFitPar() == 0; }

	Phast& _phast;

	TH1D* _hStatistics;

	// output tree variables
	int    _runNumber;
	int    _spillNumber;
	int    _eventInSpill;
	int    _triggerMask;
	int    _masterTrigger;
	int    _DT0;
	double _tcsPhase;
	double _timeInSpill;
	int    _isAcceptedMC;

};

#endif  // SELECTIONHELPER_H

#ifndef RECOILDEFINITION_H
#define RECOILDEFINITION_H

#include "Phast.h"

class RecoilDefinition
{
public:

	RecoilDefinition(TTree& tree);
	~RecoilDefinition() { }

	void fillMC(const PaEvent& event, const PaMCvertex& vertex);

	const PaMCtrack& getRecoilTrackMC(const PaEvent& event, const PaMCvertex& vertex);

private:

	Phast& _phast;

	static const int MC_PID_PROTON = 14;  // proton

	// output tree variables
	double _momentumMCX;
	double _momentumMCY;
	double _momentumMCZ;

};

#endif  // RECOILDEFINITION_H

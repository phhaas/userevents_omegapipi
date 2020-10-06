#ifndef PHOTONDEFINITION_H
#define PHOTONDEFINITION_H

#include "Phast.h"

class PhotonDefinition
{

public:

	PhotonDefinition(TTree& tree);
	~PhotonDefinition() { }

	void fill(const vector<TLorentzVector>& photons, const PaVertex& vertex) { }
	void fillMC(const PaEvent& event, const PaMCvertex& vertex);


private:

	void fillMCPairsOmegaPi0(const PaEvent& event, const PaMCvertex& vertex);

	Phast& _phast;

	static const int MAX_CLUSTERS  = 200;
	static const int MC_PID_PHOTON = 1;  // photon
	static const int MC_PID_PI0    = 7;  // Pi0
	static const int MC_PID_PROTON = 14;  // proton
	static const int MC_PID_OMEGA  = 24;  // OMEGA

	// output tree variables
	int                 _numberPhotonsMC;
	std::vector<double> _momentumMCX;
	std::vector<double> _momentumMCY;
	std::vector<double> _momentumMCZ;

};

#endif  // PHOTONDEFINITION_H

#ifndef SCATTEREDDEFINITION_H
#define SCATTEREDDEFINITION_H

#include "Phast.h"

class ScatteredDefinition
{

public:

	ScatteredDefinition(TTree& tree, const int numberChargedPartricles, const int numberNeutralPartricles);
	~ScatteredDefinition() { }

	void fill  (const PaEvent& event, const vector<PaParticle>& particles, const PaVertex& vertex);
	void fillMC(const vector<PaMCtrack>& tracks, const PaMCvertex& vertex);

	vector<PaParticle> getScatteredParticles       (const PaEvent& event, const PaVertex&   vertex);
	vector<PaTrack>    getScatteredParticleTracks  (const PaEvent& event, const PaVertex&   vertex);
	vector<PaMCtrack>  getScatteredParticleMCtracks(const PaEvent& event, const PaMCvertex& vertex);

private:

	int getNumberTotalParticles()   const { return _numberTotalParticles;   }
	int getNumberChargedParticles() const { return _numberChargedParticles; }
	int getNumberNeutralParticles() const { return _numberNeutralParticles; }

	const vector<double>& getMomentumX()   const { return _momentumX;   }
	const vector<double>& getMomentumY()   const { return _momentumY;   }
	const vector<double>& getMomentumZ()   const { return _momentumZ;   }
	const vector<double>& getMomentumMCX() const { return _momentumMCX; }
	const vector<double>& getMomentumMCY() const { return _momentumMCY; }
	const vector<double>& getMomentumMCZ() const { return _momentumMCZ; }

	Phast&  _phast;

	static const int MC_PID_PI_POS = 9;
	static const int MC_PID_PI_NEG = 8;
	static const int MC_PID_PI_0   = 7;
	static const int MC_PID_PROTON = 14;

	const int _numberChargedParticles;
	const int _numberNeutralParticles;
	int       _numberTotalParticles;
	
	// output tree variables
	std::vector<double> _momentumX;
	std::vector<double> _momentumY;
	std::vector<double> _momentumZ;
	std::vector<int>    _charge;
	std::vector<int>    _pid;
	std::vector<double> _minimalZ;
	std::vector<double> _maximalZ;
	std::vector<double> _time;
	std::vector<double> _chi2;
	std::vector<double> _XX0;
	std::vector<double> _momentumMCX;
	std::vector<double> _momentumMCY;
	std::vector<double> _momentumMCZ;
	std::vector<int>    _chargeMC;
	std::vector<int>    _pidMC;

};

#endif  // SCATTEREDDEFINITION_H

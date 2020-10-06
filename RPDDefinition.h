#ifndef RPDDEFINITIONS_H
#define RPDDEFINITIONS_H

#include "TF1.h"
#include "TH1D.h"
#include "TSpline.h"
#include "TString.h"

#include "Phast.h"

#include "RPD_Helper.h"

class RPDDefinition
{

public:

	RPDDefinition(TTree& tree);
	~RPDDefinition() { }

	void fill(const PaEvent& event, const PaVertex& vertex);

	bool cutHasTracks    (const PaEvent& event, const PaVertex& vertex) const;
	bool cutHasBestProton(const PaEvent& event, const PaVertex& vertex) const;

private:

	const vector<TLorentzVector>& getProtons   () const { return _rpd.vTracks(); }
	const TLorentzVector&         getBestProton() const;

	int    getNumberTracks()    const { return _numberTracks;    }
	bool   hasTracks()          const { return _hasTracks == 1;  }
	int    getIndexBestProton() const { return _indexBestProton; }
	const TVector3& getMomentumVector() const { return *_momentumVector; }
	double getMomentumX() const { return getMomentumVector().X(); }
	double getMomentumY() const { return getMomentumVector().Y(); }
	double getMomentumZ() const { return getMomentumVector().Z(); }
	double getBeta()  const { return _beta; }
	double getPhi()   const;
	double getTheta() const;
	double getEnergy         (const int       useMethod) const;  // 0: no correction, 1: old, 2: new
	double getCorrectedEnergy(const int       useMethod,   // 0: no correction, 1: old, 2: new; useRun: lead data run 80608 (16.10.2009 13:39), ended with run 81093 (28.10.2009 07:56)
	                          const PaEvent&  event,
	                          const PaVertex& vertex) const;

	double getHitTime()           const { return _hitTime;           }
	int    getHitsRingA()         const { return _hitsRingA;         }
	int    getHitsRingB()         const { return _hitsRingB;         }
	double getZHitPosition()      const { return _zHitPosition;      }
	double getZHitPositionRingA() const { return _zHitPositionRingA; }
	double getZHitPositionRingB() const { return _zHitPositionRingB; }
	double getEnergyLossRingA()   const { return _energyLossRingA;   }
	double getEnergyLossRingB()   const { return _energyLossRingB;   }

	double getCorrectedEnergyNew(const double pRPD, const int materialType, const double passedThickness) const;
	double getInterpolatedEnergy(const double pRPD, const int materialType, const double passedThickness) const;
	double interpolate3(const double xx1, const double xx2, const double xx3, const double yy1, const double yy2, const double yy3, const double x) const;

	Phast& _phast;

	static const bool RPD_CORRECT_MOMENTUM = false;
	static const bool RPD_CORRECT_ANGLES   = false;
	RPD& _rpd;

	static constexpr double PROTON_MASS = 0.9382720813;  // [GeV]
	std::vector<TLorentzVector> _protons;
	TLorentzVector              _bestProton;
	TVector3*                   _momentumVector;

	// output tree variables
	int    _numberTracks;
	int    _hasTracks;
	int    _indexBestProton;
	double _momentumX;
	double _momentumY;
	double _momentumZ;
	double _beta;
	double _phi;
	double _theta;
	double _energyNoCorrection;
	double _energyOldCorrection;
	double _energyNewCorrection;
	double _hitTime;
	int    _hitsRingA;
	int    _hitsRingB;
	double _zHitPosition;
	double _zHitPositionRingA;
	double _zHitPositionRingB;
	double _energyLossRingA;
	double _energyLossRingB;
	bool   _isDT0;

};

#endif  // RPDDEFINITIONS_H

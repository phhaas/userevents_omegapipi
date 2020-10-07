#ifndef BEAMDEFINITION_H
#define BEAMDEFINITION_H

#include "Phast.h"

class BeamDefinition
{

	public:

		BeamDefinition(TTree& tree);
		~BeamDefinition() { }

		void fill  (const PaEvent&    event);
		void fillMC(const PaMCtrack&  track);

		const double& getGradientdXdZ() const { return _gradX;  }
		const double& getGradientdYdZ() const { return _gradY;  }
		const double& getTime()         const { return _time;   }
		const int&    getCharge()       const { return _charge; }
		const int&    getPID()          const { return _pid;    }

		const PaParticle& getBeamParticle(const PaEvent& event, const PaVertex&   vertex) { return event.vParticle(vertex.InParticle());                        }
		const PaTrack&    getBeamTrack   (const PaEvent& event, const PaVertex&   vertex) { return event.vTrack(event.vParticle(vertex.InParticle()).iTrack()); }
		const PaMCtrack&  getBeamTrackMC (const PaEvent& event, const PaMCvertex& vertex);

	private:

		Phast& _phast;

		static const int MC_PID_BEAM = 44;  // pi- beam

		// output tree variables
		double _gradX;
		double _gradY;
		double _time;
		double _firstZ;
		double _lastZ;
		double _minimalZ;
		double _maximumZ;
		double _chi2;
		int    _charge;
		int    _pid;
		double _momentumMCX;
		double _momentumMCY;
		double _momentumMCZ;

};

#endif  // BEAMDEFINITION_H

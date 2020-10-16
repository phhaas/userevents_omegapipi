#ifndef ECAL_DEFINITION_H
#define ECAL_DEFINITION_H

#include "Phast.h"

class ECALDefinition
{

	public:

		ECALDefinition(TTree& tree);
		~ECALDefinition() { }

		void fill(const PaEvent& event);

		unsigned long getNumberNeutralClusters(const PaEvent& event) { return getVectorNeutrals(event).size(); }

	private:

		const std::vector<int>& getVectorNeutrals(const PaEvent& event);
		struct coords2D {double X,Y;};
		void getIsolation(const PaEvent& event);

		Phast& _phast;

		const PaCalorimeter& _ECAL1;
		const PaCalorimeter& _ECAL2;
		std::vector<int>     _vectorNeutrals;  // particle indices of photon candidates
		std::vector<coords2D> _chargedCoordECAL1;
		std::vector<coords2D> _chargedCoordECAL2;

		// Both ECALs
		unsigned long       _clusterNumber;
		std::vector<int> _clusterIndex;
		std::vector<double> _clusterX;
		std::vector<double> _clusterY;
		std::vector<double> _clusterZ;
		std::vector<double> _clusterXError;
		std::vector<double> _clusterYError;
		std::vector<double> _clusterZError;
		std::vector<double> _clusterEnergy;
		std::vector<double> _clusterEnergyError;
		std::vector<double> _clusterTime;
		std::vector<int> _clusterSize;
		std::vector<int>    _clusterIndexCell;
		std::vector<double> _clusterXInCell;
		std::vector<double> _clusterYInCell;
		std::vector<double> _clusterIsolation;

};

#endif  // ECAL_DEFINITION_H

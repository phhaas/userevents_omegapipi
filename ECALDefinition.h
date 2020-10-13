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
		unsigned long       _neutralClusterNumber;
		std::vector<int> _neutralClusterIndex;
		std::vector<double> _neutralClusterX;
		std::vector<double> _neutralClusterY;
		std::vector<double> _neutralClusterZ;
		std::vector<double> _neutralClusterXError;
		std::vector<double> _neutralClusterYError;
		std::vector<double> _neutralClusterZError;
		std::vector<double> _neutralClusterEnergy;
		std::vector<double> _neutralClusterEnergyError;
		std::vector<double> _neutralClusterTime;
		std::vector<int> _neutralClusterSize;
		std::vector<int>    _neutralClusterIndexCell;
		std::vector<double> _neutralClusterXInCell;
		std::vector<double> _neutralClusterYInCell;
		std::vector<double> _neutralClusterIsolation;

};

#endif  // ECAL_DEFINITION_H

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

		Phast& _phast;

		const PaCalorimeter& _ECAL1;
		const PaCalorimeter& _ECAL2;
		std::vector<int>     _vectorNeutrals;  // particle indices of photon candidates

		// output tree variables
		//TODO are separate sets for ECAL 1 and ECAL 2 really needed?
		//     wouldn't it be easier to add a vector with the ECAL indices?
		// ECAL 1
		unsigned long       _ECAL1ClusterNumber;
		// ECAL 2
		unsigned long       _ECAL2ClusterNumber;
		std::vector<int>    _ECAL2IndexCell;
		std::vector<double> _ECAL2ClusterXInCell;
		std::vector<double> _ECAL2ClusterYInCell;
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

};

#endif  // ECAL_DEFINITION_H

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

		const std::vector<double>& getECAL1clusterX() const { return _ECAL1ClusterX; };
		const std::vector<double>& getECAL1clusterY() const { return _ECAL1ClusterY; };
		const std::vector<double>& getECAL1clusterZ() const { return _ECAL1ClusterZ; };

		unsigned long getECAL1ClusterNumber(const PaEvent& event);
		unsigned long getECAL2ClusterNumber(const PaEvent& event);
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
		std::vector<double> _ECAL1ClusterX;
		std::vector<double> _ECAL1ClusterY;
		std::vector<double> _ECAL1ClusterZ;
		std::vector<double> _ECAL1ClusterXError;
		std::vector<double> _ECAL1ClusterYError;
		std::vector<double> _ECAL1ClusterZError;
		std::vector<double> _ECAL1ClusterEnergy;
		std::vector<double> _ECAL1ClusterEnergyError;
		std::vector<double> _ECAL1ClusterTime;
		std::vector<int>    _ECAL1ClusterSize;
		// ECAL 2
		unsigned long       _ECAL2ClusterNumber;
		std::vector<double> _ECAL2ClusterX;
		std::vector<double> _ECAL2ClusterY;
		std::vector<double> _ECAL2ClusterZ;
		std::vector<double> _ECAL2ClusterXError;
		std::vector<double> _ECAL2ClusterYError;
		std::vector<double> _ECAL2ClusterZError;
		std::vector<double> _ECAL2ClusterEnergy;
		std::vector<double> _ECAL2ClusterEnergyError;
		std::vector<double> _ECAL2ClusterTime;
		std::vector<int>    _ECAL2ClusterSize;
		std::vector<int>    _ECAL2IndexCell;
		std::vector<double> _ECAL2ClusterXInCell;
		std::vector<double> _ECAL2ClusterYInCell;
		// Both ECALs
		unsigned long       _neutralClusterNumber;
		std::vector<double> _neutralClusterX;
		std::vector<double> _neutralClusterY;
		std::vector<double> _neutralClusterZ;
		std::vector<double> _neutralClusterXError;
		std::vector<double> _neutralClusterYError;
		std::vector<double> _neutralClusterZError;
		std::vector<double> _neutralClusterEnergy;
		std::vector<double> _neutralClusterEnergyError;
		std::vector<double> _neutralClusterTime;
		std::vector<double> _neutralClusterSize;

};

#endif  // ECAL_DEFINITION_H
